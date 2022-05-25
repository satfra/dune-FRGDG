#install headers
install(FILES driver.hh
              backupRestoreFacility.hh
              DGTimeController.hh
              explicitonestepLDG.hh
              getIntersection.hh
              instationaryPDE.hh
              ovlpistlsolverbackend_extension.hh
              rejectStep.hh
              stationaryPDE.hh
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/FRGDG/solver)
