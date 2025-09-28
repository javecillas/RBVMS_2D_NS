cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
project(RBVMS_2D_Project)

# Eigen3 include path provided by user via GUI or -D flag
set(EIGEN3_INCLUDE_DIR "" CACHE PATH "Path to Eigen3 include folder")

if(EIGEN3_INCLUDE_DIR)
  message(STATUS "Using user-provided Eigen3 path: ${EIGEN3_INCLUDE_DIR}")
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message(STATUS "No Eigen3 path provided. You can set EIGEN3_INCLUDE_DIR in CMake GUI or with -D flag.")
endif()

# Find the Boost library
find_package(Boost REQUIRED QUIET COMPONENTS system thread)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
else()
  message(FATAL_ERROR "Boost not found! ${INSTALL_MESSAGE}")
endif()

# Specify where the headers files must be found
include_directories(
  "${PROJECT_SOURCE_DIR}/analysis"
  "${PROJECT_SOURCE_DIR}/error"
  "${PROJECT_SOURCE_DIR}/material"
  "${PROJECT_SOURCE_DIR}/mesh"
  "${PROJECT_SOURCE_DIR}/pre"
)

# Define where cmake should look for other CMakeList.txt files
add_subdirectory(analysis)
add_subdirectory(error)
add_subdirectory(material)
add_subdirectory(mesh)
add_subdirectory(pre)

# Here is where we define and create the library
add_subdirectory(lib)

# Add the executable
add_executable(rbvms2d main.cpp)

# Link the libraries
target_link_libraries(rbvms2d fem_lib ${Boost_LIBRARIES})