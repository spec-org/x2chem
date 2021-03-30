add_library(x2chem STATIC x2chem.cpp)
target_include_directories(x2chem PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(x2chem PUBLIC X2Chem::blaspp X2Chem::lapackpp)
