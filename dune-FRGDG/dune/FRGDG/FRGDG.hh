#pragma once

#include <exception>
#include <string>

// DUNE includes:
// common module
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/exceptions.hh>

// current module
#include <dune/FRGDG/simulationmgr.hh>
namespace Dune
{
template<class SIM>
  void startFRGSimulation(int argc, char **argv, std::string iniName)
  {
    try
    {
      // we open the configuration file through the parser and store the results in a ParameterTree
      ParameterTreeParser ptreeparser;
      ParameterTree ptree;
      ptreeparser.readINITree(iniName, ptree);
      ptreeparser.readOptions(argc, argv, ptree);

      // Initialize MPI for use with threaded applications - this is a custom Dune extension for dune_FRGDG
      const MPIHelper &helper = Dune::MPIHelper::instance_threaded(argc, argv);

      // If the program is running on a singular node, just start the simulation manager on this node
      if (MPIHelper::isFake)
      {
        std::cout << "This is a single-node program.\n";
        SimulationMgr<SIM> mgr(ptree,helper);
        mgr.run();
      }
      // If the program is running on a multiple nodes, but we do only a single simulation, run the simulation manager on all nodes
        else
        {
          if(helper.rank() == 0)
            std::cout << "This is a cluster distributed program, running on " << helper.size() << " nodes.\n";
          SimulationMgr<SIM> mgr(ptree,helper);
          mgr.run();
        }
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
}
