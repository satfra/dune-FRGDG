FILE(GLOB SRCFILES ./common/*.cc )
set(ExPath "executables")
set(IniPath "ini")
find_package(GSL REQUIRED)

macro(generate_executable name)
  execute_process(COMMAND bash ./generateSource.sh ${name} ${ExPath} ${IniPath}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${IniPath}/${name}.ini" "${CMAKE_CURRENT_BINARY_DIR}/${name}.ini")
  add_executable("${name}" ${ExPath}/${name}.cc ${SRCFILES})
  target_link_dune_default_libraries("${name}")
  target_link_libraries("${name}" ${GSL_LIBRARIES})
endmacro()

generate_executable(anharmonicOscillator)
generate_executable(largeN)
generate_executable(ON)
