cmake_minimum_required(VERSION 3.7)
project(kalman_filters LANGUAGES CXX)
get_filename_component(SOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, TPCSOLUTION_DIR: ${SOLUTION_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS
    include/kalman_filter_debug.h
    include/kalman_filter_linear.h

    include/kalman_filter_extended.h
    include/kalman_filter_extended_square_root.h

    include/kalman_filter_cubature.h
    include/kalman_filter_cubature_square_root.h
    include/kalman_filter_extended_cubature.h
    include/kalman_filter_extended_cubature_square_root.h

    include/kalman_filter_unscented.h
    include/kalman_filter_unscented_square_root.h
    include/kalman_filter_extended_unscented.h
    include/kalman_filter_extended_unscented_square_root.h

    ../3d_party/timing.h
    ../3d_party/qr_decomposition.h
    ../3d_party/blas0.h
    ../3d_party/blas1_d.h
    ../3d_party/linpack_d.h
    ../3d_party/cholupdate_linpack.h
)

set(SOURCES
    src/kalman_filter_debug.cpp

    ../3d_party/timing.cpp
    ../3d_party/qr_decomposition.cpp
    ../3d_party/blas0.cpp
    ../3d_party/blas1_d.cpp
    ../3d_party/linpack_d.cpp
    ../3d_party/cholupdate_linpack.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES}) # Статическая библиотека

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/../3d_party
)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}
    PUBLIC
        ${ARMADILLO_LIBRARIES}
)

# Включить замеры времени сборки
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CMAKE_COMMAND} -E time")
