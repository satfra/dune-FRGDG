#install headers
install(FILES modelInterface.hh
              modelInterfaceLDG.hh
              anharmonicOscillator.hh
              largeN.hh
              ON.hh
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/FRGDG/model)

