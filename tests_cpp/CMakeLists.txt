add_executable(velesquant_tests
    test_sabr.cpp
    test_hw.cpp
    test_tree.cpp
    test_pde.cpp
    test_utilities.cpp
)

target_link_libraries(velesquant_tests
    PRIVATE
    GTest::gtest_main
    velesquant_core
)

include(GoogleTest)
gtest_discover_tests(velesquant_tests)
