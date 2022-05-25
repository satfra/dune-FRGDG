#install headers
install(FILES fluxInterface.hh
              centralFlux.hh
              LLFFlux.hh
              upwindFlux.hh
              godunovFlux.hh
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/FRGDG/numericalflux)

