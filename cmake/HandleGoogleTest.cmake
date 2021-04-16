include ( FetchContent )

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
)

FetchContent_MakeAvailable( googletest )

add_library( X2Chem::gtest ALIAS gtest )
add_library( X2Chem::gmock ALIAS gmock )
