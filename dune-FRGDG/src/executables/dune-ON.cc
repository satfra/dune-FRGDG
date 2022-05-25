#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/FRGDG/FRGDG.hh>
#include <dune/FRGDG/simulation.hh>
#include <dune/FRGDG/model/ON.hh>

int main(int argc, char **argv)
{
  using SIM = Simulation<ON::SimSet>;
  Dune::startFRGSimulation<SIM>(argc, argv, "ON.ini");
  return 0;
}
