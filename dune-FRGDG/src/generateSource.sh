#!/bin/bash

name=$1
ExPath=$2
IniPath=$3

echo ${ExPath}

mkdir -p ${ExPath}

echo "#ifdef HAVE_CONFIG_H
#include \"config.h\"
#endif

#include <dune/FRGDG/FRGDG.hh>
#include <dune/FRGDG/simulation.hh>
#include <dune/FRGDG/model/${name}.hh>

int main(int argc, char **argv)
{
  using SIM = Simulation<${name}::SimSet>;
  Dune::startFRGSimulation<SIM>(argc, argv, \"${name}.ini\");
  return 0;
}" > ${ExPath}/${name}.cc
touch ${IniPath}/${name}.ini
