file(GLOB TEST_FILES test*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)


add_executable(unit_tests catch_main.cpp ${TEST_FILES})

# Properly link your library target
target_link_libraries(unit_tests PRIVATE HydroGravLib)

# add tests

add_test(NAME UnitTests COMMAND unit_tests)