cmake_minimum_required(VERSION 3.7)
project(kalman_filters LANGUAGES CXX)
get_filename_component(SOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, TPCSOLUTION_DIR: ${SOLUTION_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(HEADERS
    kalman_debug.h
    kalman_base.h
    # kalman_lkf.h

    kalman_ekf.h
    kalman_srekf.h

    kalman_ckf.h
    kalman_srckf.h
    kalman_eckf.h
    kalman_sreckf.h

    kalman_ukf.h
    kalman_srukf.h
    kalman_eukf.h
    kalman_sreukf.h
    
    ../3d_party/qr_decomposition.h    
)

set(SOURCES
    ../3d_party/qr_decomposition.cpp    
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES}) # Статическая библиотека

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
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
