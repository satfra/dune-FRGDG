#install headers
install(FILES FRGDG.hh
              simulation.hh
              simulationInterface.hh
              simulationmgr.hh
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/FRGDG)

add_subdirectory(common)
add_subdirectory(solver)
add_subdirectory(numericalflux)
add_subdirectory(model)
add_subdirectory(localoperator)
