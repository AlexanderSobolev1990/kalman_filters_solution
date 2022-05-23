cmake_minimum_required(VERSION 3.7)
project(compare LANGUAGES CXX)
get_filename_component(SOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, SOLUTION_DIR: ${SOLUTION_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(BLAS REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS system program_options)
include_directories(/usr/include/python3.8)

add_executable(${PROJECT_NAME}
    compare.h
    compare.cpp
    main.cpp
    )
target_link_libraries(${PROJECT_NAME}
    ${ARMADILLO_LIBRARIES}
    kalman_filters
    libpython3.8.so
    Boost::system
    Boost::program_options
    )
target_include_directories(${PROJECT_NAME}
    PUBLIC        
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${SOLUTION_DIR}/../kalman_filters/include
    )
