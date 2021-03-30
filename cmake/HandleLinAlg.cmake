include ( FetchContent )

FetchContent_Declare(
  blaspp
  GIT_REPOSITORY https://bitbucket.org/icl/blaspp.git
)

FetchContent_Declare(
  lapackpp 
  GIT_REPOSITORY https://bitbucket.org/icl/lapackpp.git
)

FetchContent_MakeAvailable( blaspp )
FetchContent_MakeAvailable( lapackpp )

add_library( X2Chem::blaspp ALIAS blaspp )
add_library( X2Chem::lapackpp ALIAS lapackpp )
