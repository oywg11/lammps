if(BUILD_TOOLS)
  add_executable(binary2txt ${LAMMPS_TOOLS_DIR}/binary2txt.cpp)
  target_compile_definitions(binary2txt PRIVATE -DLAMMPS_${LAMMPS_SIZES})
  install(TARGETS binary2txt DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(stl_bin2txt ${LAMMPS_TOOLS_DIR}/stl_bin2txt.cpp)
  install(TARGETS stl_bin2txt DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(reformat-json ${LAMMPS_TOOLS_DIR}/json/reformat-json.cpp)
  target_include_directories(reformat-json PRIVATE ${LAMMPS_SOURCE_DIR})
  install(TARGETS reformat-json DESTINATION ${CMAKE_INSTALL_BINDIR})

  include(CheckGeneratorSupport)
  if(CMAKE_GENERATOR_SUPPORT_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
      enable_language(Fortran)
      add_executable(chain.x ${LAMMPS_TOOLS_DIR}/chain.f90)
      target_link_libraries(chain.x PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
      add_executable(micelle2d.x ${LAMMPS_TOOLS_DIR}/micelle2d.f90)
      target_link_libraries(micelle2d.x PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
      install(TARGETS chain.x micelle2d.x DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
      message(WARNING "No suitable Fortran compiler found, skipping build of 'chain.x' and 'micelle2d.x'")
    endif()
  else()
    message(WARNING "CMake build doesn't support Fortran, skipping build of 'chain.x' and 'micelle2d.x'")
  endif()

  enable_language(C)
  get_filename_component(MSI2LMP_SOURCE_DIR ${LAMMPS_TOOLS_DIR}/msi2lmp/src ABSOLUTE)
  file(GLOB MSI2LMP_SOURCES CONFIGURE_DEPENDS ${MSI2LMP_SOURCE_DIR}/[^.]*.c)
  add_executable(msi2lmp ${MSI2LMP_SOURCES})
  if(STANDARD_MATH_LIB)
    target_link_libraries(msi2lmp PRIVATE ${STANDARD_MATH_LIB})
  endif()
  install(TARGETS msi2lmp DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${LAMMPS_DOC_DIR}/msi2lmp.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

  add_subdirectory(${LAMMPS_TOOLS_DIR}/phonon ${CMAKE_BINARY_DIR}/phana_build)
endif()

if(BUILD_LAMMPS_GUI)
  include(ExternalProject)
  option(LAMMPS_GUI_USE_PLUGIN "Load LAMMPS library dynamically at runtime" OFF)
  option(LAMMPS_GUI_BUILD_DOCS "Build LAMMPS-GUI HTML documentation" OFF)
  ExternalProject_Add(lammps-gui_build
    GIT_REPOSITORY https://github.com/akohlmey/lammps-gui.git
    GIT_TAG main
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    CMAKE_ARGS -D BUILD_DOC=${LAMMPS_GUI_BUILD_DOCS}
               -D LAMMPS_GUI_USE_PLUGIN=${LAMMPS_GUI_USE_PLUGIN}
               -D LAMMPS_SOURCE_DIR=${LAMMPS_SOURCE_DIR}
               -D LAMMPS_LIBRARY=$<TARGET_FILE:lammps>
               -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -D CMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -D CMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
               -D CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    DEPENDS lammps
    BUILD_BYPRODUCTS <INSTALL_DIR>/bin/lammps-gui
  )
endif()
