#include <gtest/gtest.h>

TEST( DUMMY_TEST, DUMMY ) {
  ASSERT_TRUE(true);
}


// Required
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
