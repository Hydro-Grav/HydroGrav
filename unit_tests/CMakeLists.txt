file(GLOB TEST_FILES test*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_executable(unit_tests catch_main.cpp ${TEST_FILES})

target_link_libraries(unit_tests PRIVATE HydroGravLib)

enable_testing()
add_test(NAME UnitTests COMMAND unit_tests)
