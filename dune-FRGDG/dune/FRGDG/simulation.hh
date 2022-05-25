#pragma once


#include <thread>
#include <future>
#include <chrono>
#include <cassert>
#include <sstream>
#include <memory>
#include <exception>

#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/common/structuredgridfactory.hh>

#include <dune/common/parametertree.hh>
#include <dune/common/parallel/mpicommunication.hh>

// this module
#include <dune/FRGDG/common/utils.hh>
#include <dune/FRGDG/common/logger.hh>

#include <dune/FRGDG/solver/driver.hh>

template<class SET>
class Simulation
{
  public:
    static constexpr int dim = SET::Traits::dim;
    using GridConstructor = typename SET::Traits::GridConstructor;
    using Grid = typename GridConstructor::Grid;
    using GV = typename GridConstructor::GV;
    using RF = typename SET::Traits::RF;
    using ResultType = Dune::FieldVector<RF,dim>;
    using CommType = Dune::Communication<MPI_Comm>;

    template<int order>
    using DRV = Driver<SET,order>;

    Simulation(Dune::ParameterTree ptree_, const CommType& mpicomm_, const unsigned group, Logger& log_) :
      result(0.), mpicomm(mpicomm_), log(log_), ptree(ptree_), running(false) {}

    const Dune::FieldVector<RF,dim>& getResult() const
    {
      return result;
    }

    void start();
    void finish();
    bool finished();

  protected:
    template<int curOrderIdx>
      void runDriver(const GV& gv, const int& givenOrder) {
        //Execute some code
        constexpr int N = SET::Traits::orders[curOrderIdx];
        if(N == givenOrder)
        {
          typename SET::Traits::template FEM<N> fem;
          DRV<N> driver(gv, fem, ptree, mpicomm,log);
          try {
            driver.Run();
          }
          catch (const Dune::Exception &e)
          {
            std::cerr << "Dune reported error: " << e << std::endl;
          }
          catch (const std::exception &exc)
          {
            std::cerr << "Exception thrown: " << exc.what() << std::endl;
          }
          result = driver.GetResult().back();
        }
        else
        {
          if constexpr (curOrderIdx > 0)
            runDriver<curOrderIdx-1>(gv, givenOrder);
          else
            throw(std::runtime_error("The given fem.degree has not been compiled into this runtime!"));
        }
      }

    void program();
    ResultType result;
    CommType mpicomm;
    Logger& log;
    Dune::ParameterTree ptree;

    std::future<void> future;
    bool running;
};

  template<class SET>
void Simulation<SET>::program()
{
  try 
  {
    if(mpicomm.rank() == 0)
      std::cout << "Simulation started\n";
    log.log("Started Simulation");

    // We construct a rectangular grid from data in the ptree
    Dune::GridPtr<Grid> grid = GridConstructor::getGrid(ptree,mpicomm);
    GV gv = grid->leafGridView();
    log.log("Created Grid");
      
    // Depending on the chosen degree of the FEM we need to specify a different Element Map
    const int degree = ptree.get<int>("fem.degree");
    log.log("Starting Driver");
    runDriver<SET::Traits::orders.size()-1>(gv, degree);
    log.log("Driver finished");

    if(mpicomm.rank() == 0)
      std::cout << "Simulation finished\n";

  } 
  catch (const Dune::Exception &e) 
  { 
    std::cerr << "Dune reported error: " << e << std::endl; 
  } 
  catch (const std::exception &exc) 
  { 
    std::cerr << "Exception thrown: " << exc.what() << std::endl; 
  } 
}

  template<class SET>
void Simulation<SET>::start()
{
  if(!running)
  {
    future = std::async(std::launch::async, &Simulation<SET>::program, this);
    running = true;
  }
}

  template<class SET>
void Simulation<SET>::finish()
{
  if(future.valid())
    future.wait();
  running = false;
}

  template<class SET>
bool Simulation<SET>::finished()
{
  if(future.valid())
    return future.wait_for(std::chrono::nanoseconds(1)) == std::future_status::ready;
  return false;
}
