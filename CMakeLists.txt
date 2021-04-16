cmake_minimum_required(VERSION 3.11)

#------------------------------------------------------------------------------
# Project declaration
#------------------------------------------------------------------------------

project(
  X2Chem
  VERSION 0.1.0
  LANGUAGES C CXX
  DESCRIPTION
  "X2Chem is a lightweight library to provide the functionality required for \
relativistic electronic structure calculations to existing quantum chemistry \
code bases."
)


# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

option(ENABLE_TESTS "Build unit tests" ON)


# -----------------------------------------------------------------------------
# Global properties
# -----------------------------------------------------------------------------

# C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)


#------------------------------------------------------------------------------
# External dependencies
#------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# SMP parallelism
find_package(OpenMP)

# Linear algebra
include(HandleLinAlg)


#------------------------------------------------------------------------------
# The actual library
#------------------------------------------------------------------------------

add_subdirectory(src)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------

if(ENABLE_TESTS)
  include(HandleGoogleTest)
  enable_testing()
  add_subdirectory(test)
endif()
