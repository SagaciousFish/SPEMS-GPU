# Add the Python script to the project
message(STATUS "Adding Python-based tests")
message(STATUS "Python executable: ${PYTHON_EXECUTABLE}")
message(STATUS "CMake Source Directory: ${CMAKE_CURRENT_SOURCE_DIR}")

set(PYTHON_TEST_SCRIPT "${CMAKE_BINARY_DIR}/test/run_profiles.py")

# Ensure Python is available
find_package(Python COMPONENTS Interpreter REQUIRED)

# Add a custom target to run the Python tests
add_custom_target(run_python_tests
        COMMAND ${PYTHON_EXECUTABLE} ${PYTHON_TEST_SCRIPT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running Python-based tests"
)

set(TEST_TARGET parEMS)

# Ensure the CLI binary is built before running the tests
add_dependencies(run_python_tests ${TEST_TARGET})

set(CLI_BINARY_PATH ${CMAKE_BINARY_DIR}/test)

configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/run_profiles.py
        ${CMAKE_BINARY_DIR}/test/run_profiles.py
        @ONLY
)

add_test(
        NAME PythonTests
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/test/run_profiles.py
)