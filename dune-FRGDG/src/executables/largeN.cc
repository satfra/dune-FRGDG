#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/FRGDG/FRGDG.hh>
#include <dune/FRGDG/simulation.hh>
#include <dune/FRGDG/model/largeN.hh>

int main(int argc, char **argv)
{
  using SIM = Simulation<largeN::SimSet>;
  Dune::startFRGSimulation<SIM>(argc, argv, "largeN.ini");
  return 0;
}
