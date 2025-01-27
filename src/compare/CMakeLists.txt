cmake_minimum_required(VERSION 3.7)
project(compare LANGUAGES CXX)
get_filename_component(SOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, SOLUTION_DIR: ${SOLUTION_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(Armadillo REQUIRED)
#find_package(Boost REQUIRED COMPONENTS system program_options)
include_directories(/usr/include/python3.8)

add_executable(${PROJECT_NAME}
    main.cpp
    compare_run.h
    compare_run_main.cpp
    compare_run_rmse.cpp        
    # 
    ../3d_party/compare.h
    ../3d_party/consts.h
    ../3d_party/convert.h
    ../3d_party/convert.cpp
    ../3d_party/geodesy.h
    ../3d_party/geodesy.cpp
    ../3d_party/units.h
    ../3d_party/radars.h
    ../3d_party/radars.cpp    
    ../3d_party/numericalmethods.h
    ../3d_party/statistics.h
    ../3d_party/statistics.cpp
    ../3d_party/timing.h
    ../3d_party/timing.cpp
    )
target_link_libraries(${PROJECT_NAME}
    ${OpenMP_CXX_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    kalman_filters
    libpython3.8.so
    /lib/x86_64-linux-gnu/libboost_system.so.1.67.0
    /lib/x86_64-linux-gnu/libboost_program_options.so.1.67.0
    )
target_include_directories(${PROJECT_NAME}
    PUBLIC        
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../3d_party
    )
