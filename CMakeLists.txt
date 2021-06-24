CMake_minimum_required(VERSION 3.1)
project(asap VERSION 1.0.0)
find_package(PythonLibs REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O2 -Wall")

add_library(asaplib
 	src/asap.cpp)

target_include_directories(asaplib PUBLIC "${PROJECT_SOURCE_DIR}/ext/nanoflann")
target_include_directories(asaplib PRIVATE ext/pybind11/include)

message("Remember to insert the dimension of data you want to study before compiling the code using CMake!")
target_compile_definitions(asaplib PUBLIC DIMENSION=3)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11")
pybind11_add_module(asap src/asap.cpp)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(asap PRIVATE asaplib ${PYTHON_LIBRARIES})

set_target_properties( asap
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python"
)