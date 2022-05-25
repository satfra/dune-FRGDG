#pragma once

#include <sys/stat.h>
#include <exception>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/solver/newton.hh>
#include <dune/pdelab/constraints/p0.hh>

#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/common/typeFactory.hh>
#include <dune/FRGDG/common/logger.hh>

#include <dune/FRGDG/solver/DGTimeController.hh>
#include <dune/FRGDG/solver/getIntersection.hh>
#include <dune/FRGDG/solver/rejectStep.hh>
#include <dune/FRGDG/solver/instationaryPDE.hh>

template <typename SET, int order>
class Driver
{ 
  private:
    static constexpr double EPSILON = 1e-10;
  public:
    using RF = typename SET::Traits::RF;
    static constexpr unsigned dim = SET::Traits::dim;
    using GV = typename SET::Traits::GridConstructor::GV;
    using FEMDG = typename SET::Traits::template FEM<order>;

    using CommType = Dune::Communication<MPI_Comm>;

    using Scheme = typename SET::Traits::template Scheme<order>;

    using INTERSEC = Dune::PDELab::IntersectionFinder<typename Scheme::DataStorage>;

  private:
    // ptree stuff
    int verbosity;

    Dune::FieldVector<RF,dim> val_intersection;
    std::vector<Dune::FieldVector<RF,dim>> results;

    std::string simulationName;
    std::string currentCalculationName;
    std::string outputFile;
    RF maxTime;
    int timesteppingOrder;
    RF timeGridDistance;
    RF alpha;
    int subsampling;
    RF minTimeStep;
    RF orig_maxTimeStep;
    RF cur_maxTimeStep;

    RF cur_prefactor;
    RF orig_prefactor;

    // non-config
    const GV& gv;
    const FEMDG& femdg;

    Dune::ParameterTree& ptree;
    CommType mpicomm;
    Logger& log;

    Scheme scheme;

  public:
    Driver(const GV& gv_, const FEMDG &femdg_, Dune::ParameterTree &_ptree, CommType mpicomm_ = CommType(), Logger& log_ = Logger::discard())
      : gv(gv_), femdg(femdg_), ptree(_ptree), mpicomm(mpicomm_), log(log_), scheme(gv, femdg, ptree, mpicomm, log)
    {
      // read and set all parameters for the simulation
      subsampling = ptree.get("output.subsampling", order);
      if(subsampling < 1)
        subsampling = order;
      simulationName = ptree.get("output.name", "output");
      currentCalculationName = simulationName + SET::Traits::template Model<0>::getName(ptree);
      outputFile = currentCalculationName + "-data.txt";

      timesteppingOrder = ptree.get("fem.torder", int(1));
      maxTime = ptree.get("time.maxTime", (RF)1.0);
      timeGridDistance = ptree.get("time.timeGrid", RF(maxTime));
      orig_prefactor = ptree.get("time.safetyfactor", RF(1.));
      cur_prefactor = orig_prefactor;

      alpha = 0;

      orig_maxTimeStep = ptree.get("time.maxTimeStep", RF(1.0));
      cur_maxTimeStep = orig_maxTimeStep;
      minTimeStep = ptree.get("time.minTimeStep", RF(1e-14));
      verbosity = ptree.get<int>("solver.verbosity", int(0));

      val_intersection = SET::Traits::template Model<0>::getEOMIntersection(ptree);
    }

    std::vector<Dune::FieldVector<RF,dim>>& GetResult()
    {
      return results;
    }

    void Run()
    {
      timeIntegration();
    }

  private:
    std::string constructName(double T, double mu, Dune::ParameterTree ptree)
    {
      std::string simulationName = ptree.get("output.name", "output");
      return simulationName + "_mu=" + std::to_string(mu) + "_T=" + std::to_string(T);
    }
    
    /*
     * Method that integrates the time evolution using a specified solver
     *
     */
    void timeIntegration()
    {
      auto vtkSequenceWriter = scheme.constructVTKWriter(currentCalculationName, subsampling);
      vtkSequenceWriter.write(0.0, Dune::VTK::ascii);

      // initialize simulation parameters
      RF time = 0.0;
      RF dt = 0.;

      // add the result at time, setting it to 0
      results.resize(1);
      results[0] = 0.;

      INTERSEC intersection(scheme.getDataStorage(), mpicomm);
      
      RF dt_old = 0.;

      const auto t0 = std::chrono::high_resolution_clock::now();
      while (time < maxTime)
      {
        const auto t1 = std::chrono::high_resolution_clock::now();

        dt = scheme.step(time, cur_prefactor);
        
        // output the current timestep statistics
        if(verbosity == -1 && mpicomm.rank() == 0)
        {
          const auto t2 = std::chrono::high_resolution_clock::now();
          const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
          const auto gduration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t0).count();
          std::cout << std::setw(12) << std::setprecision(4) << std::scientific 
            << "from t = " << time << " | dt = " << dt << " | to t = " << time+dt << " | calc_dt = " 
            << duration << "ms | calc_t = " << utils::timeFormat(gduration) << std::endl;
        }

        // whenever we are at a time grid point, save the current data
        if(scheme.isAtTimeGridPoint())
        {
          if(scheme.checkConsistency())
          {
            time += dt;
            if(scheme.isAtTimeGridPoint())
              saveData(vtkSequenceWriter, time, intersection, dt_old); 
            scheme.accept();
            dt_old = dt;
          }
          else 
          {
            std::cerr << "Failure. Aborting integration\n";
            saveData(vtkSequenceWriter, time, intersection, dt_old); 
            break;
          }
        }
        else
        {
          time += dt;
          scheme.accept();
          dt_old = dt;
        }
        
        // respect the break condition when the time step gets too small
        if(dt < minTimeStep)
        {
          if(verbosity == -1 && mpicomm.rank() == 0)
            std::cout << "Timestep got smaller than " << minTimeStep << ", ending the simulation\n";
          break;
        }
      }
      // save the very last data point, even if it is not a grid point
      saveData(vtkSequenceWriter, time, intersection, dt); 
    }

    /*
     * A helper function to save Data to file and attach results to internal cache
     */
    void saveData(Dune::VTKSequenceWriter<GV>& vtkSequenceWriter, RF time, INTERSEC& intersection, RF dt)
    {
      // write the solution to a vtk file
      vtkSequenceWriter.write(time, Dune::VTK::ascii);
      // get the current intersection value
      auto r = intersection.getIntersection(val_intersection);
      // and store the result vector
      results.push_back(r);
      alpha = scheme.getAlpha();
      RF max_alpha = mpicomm.max(alpha);
      
      // if write the overall information to an output file (only one node in the mpi group)
      if(mpicomm.rank() == 0)
      {
        std::ofstream outputFileStream;
        outputFileStream.open(outputFile, std::ios_base::app);
        outputFileStream << time << " " << max_alpha;
        for(unsigned i = 0; i < dim; ++i)
          outputFileStream << " " << r[i];
        outputFileStream << " " << dt << "\n";
        outputFileStream.close();
      }
    }

    /*
     * A helper function to delete all cached result data after a certain time point
     * @param time after which one should delete all cached results
     */
    void removeAfterTime(RF time)
    {
      int idx = int(std::floor(time/timeGridDistance))+1;
      results.erase(results.begin()+idx, results.end());
    }
};
