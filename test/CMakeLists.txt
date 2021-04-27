function(add_x2chem_test testname filter)
  add_test(NAME ${testname} COMMAND x2ctest --gtest_filter=${filter})
endfunction()

#------------------------------------------------------------------------------
# Basic testing executable
#------------------------------------------------------------------------------

add_executable(x2ctest x2ctest.cxx ortho.cxx)
target_link_libraries(x2ctest PUBLIC x2chem X2Chem::gtest)

add_x2chem_test("Orthonormalize" "Orthonormal.*")
