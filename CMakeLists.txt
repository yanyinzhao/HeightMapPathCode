# Template Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( HeightMapPathCode )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL REQUIRED )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )
find_package(Boost REQUIRED COMPONENTS system program_options)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

project(main C CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_executable(main src/main.cpp)
target_link_libraries(main Boost::system Boost::program_options)
target_link_libraries(main CGAL::CGAL)


