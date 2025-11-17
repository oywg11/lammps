cmake_minimum_required(VERSION 3.20)
project(testlmp CXX)

find_package(LAMMPS REQUIRED)
add_executable(testlmp main.cpp)
if(USE_MPI)
  target_compile_definitions(testlmp PRIVATE USE_MPI)
endif()
target_link_libraries(testlmp PRIVATE LAMMPS::LAMMPS)
