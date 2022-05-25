#install headers
install(FILES DG_LOP_ConEq.hh
              DG_TLOP_ConEq.hh
              LDG_LOP_ConEq.hh
              LDG_LOP_ConEq_stat.hh
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/FRGDG/localoperator)
