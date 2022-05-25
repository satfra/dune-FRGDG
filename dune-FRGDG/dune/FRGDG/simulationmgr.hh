#pragma once

// std includes
#include <vector>
#include <thread>
#include <chrono>
#include <memory>
#include <sys/stat.h>

// DUNE includes:
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/common/backuprestore.hh>

#include <dune/common/parametertree.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/mcmgmapper.hh>

// this module
#include <dune/FRGDG/common/logger.hh>
#include <dune/FRGDG/common/utils.hh>

template<class SIM>
class SimulationMgr
{
  public:
    //using GridType = Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>;
    using RF = typename SIM::RF;

    SimulationMgr(Dune::ParameterTree ptree_, const Dune::MPIHelper& mpihelper_);
    static constexpr unsigned dim = SIM::dim;

    void run();

  private:
    Dune::ParameterTree ptree;

    std::string name;

    const Dune::MPIHelper& mpihelper;
};

  template<class SIM>
SimulationMgr<SIM>::SimulationMgr(Dune::ParameterTree ptree_, const Dune::MPIHelper& mpihelper_) : 
  ptree(ptree_), mpihelper(mpihelper_) {}

std::string constructPath(Dune::ParameterTree& ptree)
{
  std::string simulationName = ptree.get("output.name", "output");
  return utils::makePath(utils::getCurrentWorkingDir()) + simulationName + "_log/";
}

template<typename Comm>
std::string constructName(Dune::ParameterTree& ptree, const Comm& comm, const unsigned group)
{
  std::string simulationName = ptree.get("output.name", "output");
  return simulationName + "_group_" + std::to_string(group) + "_node_" + std::to_string(comm.rank());
}

  template<class SIM>
void SimulationMgr<SIM>::run()
{
  Logger log(constructPath(ptree), constructName(ptree, mpihelper.getCommunication(), 0));
    
  SIM sim(ptree, mpihelper.getCommunication(), 0, log);

  const auto t1 = std::chrono::high_resolution_clock::now();
  sim.start();
  sim.finish();
  const auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

  if(mpihelper.rank() == 0)
    std::cout << "\nEvaluation time : " << duration << "s\n---------------------------------------------------\n";
}

