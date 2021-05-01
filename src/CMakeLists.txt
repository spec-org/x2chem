#
# Main library
#

add_library(x2chem STATIC x2chem.cpp)
target_link_libraries(x2chem PUBLIC X2Chem::blaspp X2Chem::lapackpp)

target_include_directories(x2chem
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

#
# Installation
#

# Library
install(TARGETS x2chem EXPORT x2chem-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
# Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Configuration files
set( confdir ${CMAKE_INSTALL_LIBDIR}/cmake/X2Chem )

install(
  EXPORT x2chem-targets
  NAMESPACE X2Chem::
  FILE X2ChemTargets.cmake
  DESTINATION ${confdir}
)

# Also export to build directory for FetchContent
export(
  EXPORT x2chem-targets
  NAMESPACE X2Chem::
  FILE ${PROJECT_BINARY_DIR}/X2ChemTargets.cmake
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/X2ChemConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/X2ChemConfig.cmake
  INSTALL_DESTINATION ${confdir}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/X2ChemConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/X2ChemConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/X2ChemConfigVersion.cmake
  DESTINATION ${confdir}
)
