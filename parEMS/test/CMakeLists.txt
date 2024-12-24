# We adding a lot of tests here
add_executable(parEMS parEMS.cpp)
target_link_libraries(parEMS ${PROJECT_NAME})

# Add the Python script to the project
add_subdirectory(python-test)
