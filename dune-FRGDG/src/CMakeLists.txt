add_subdirectory(ini)

FILE(GLOB SRCFILES ./common/*.cc )

set(ExPath "executables")
set(IniPath "ini")

add_executable("dune-anharmonicOscillator" ${ExPath}/dune-anharmonicOscillator.cc ${SRCFILES})
target_link_dune_default_libraries("dune-anharmonicOscillator")

add_executable("dune-largeN" ${ExPath}/dune-largeN.cc ${SRCFILES})
target_link_dune_default_libraries("dune-largeN")

add_executable("dune-ON" ${ExPath}/dune-ON.cc ${SRCFILES})
target_link_dune_default_libraries("dune-ON")
