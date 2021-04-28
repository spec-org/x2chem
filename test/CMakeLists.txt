function(add_x2chem_test testname filter)
  add_test(NAME ${testname} COMMAND x2ctest --gtest_filter=${filter})
endfunction()

#------------------------------------------------------------------------------
# Basic testing executable
#------------------------------------------------------------------------------

add_executable(x2ctest x2ctest.cxx boettger.cxx form1esoc.cxx ortho.cxx x2cham.cxx)
target_link_libraries(x2ctest PUBLIC x2chem X2Chem::gtest)

add_x2chem_test("Boettger" "Boettger.*")
add_x2chem_test("Form1eSOC" "Form1eSOC.*")
add_x2chem_test("Orthonormalize" "Orthonormal.*")
add_x2chem_test("X2C_Hamiltonian" "X2C_Hamiltonian.*")
