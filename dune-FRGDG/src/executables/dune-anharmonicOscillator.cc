#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/FRGDG/FRGDG.hh>
#include <dune/FRGDG/simulation.hh>
#include <dune/FRGDG/model/anharmonicOscillator.hh>

int main(int argc, char **argv)
{
  using SIM = Simulation<anharmonicOscillator::SimSet>;
  Dune::startFRGSimulation<SIM>(argc, argv, "anharmonicOscillator.ini");
  return 0;
}
