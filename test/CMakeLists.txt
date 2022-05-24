cmake_minimum_required(VERSION 3.7)
project(test_kalman_filters LANGUAGES CXX)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}")

enable_testing(true)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Boost COMPONENTS unit_test_framework REQUIRED)
find_package(BLAS REQUIRED)
find_package(Armadillo REQUIRED)
include_directories(/usr/include/python3.8)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#-----------------------------------------------------------------------------------------------------------------------
# test_kalman_filters
add_executable(test_kalman_filters test_kalman_filters.cpp)
add_test(NAME test_kalman_filters COMMAND test_kalman_filters)
target_link_libraries(test_kalman_filters kalman_filters libpython3.8.so ${Boost_LIBRARIES})
target_include_directories(test_kalman_filters PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/3d_party
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/compare
    )
#-----------------------------------------------------------------------------------------------------------------------
# test_qr_decompositions
add_executable(test_qr_decompositions test_qr_decompositions.cpp)
add_test(NAME test_qr_decompositions COMMAND test_qr_decompositions)
target_link_libraries(test_qr_decompositions kalman_filters libpython3.8.so ${Boost_LIBRARIES})
target_include_directories(test_qr_decompositions PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/3d_party
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/compare
    )
#-----------------------------------------------------------------------------------------------------------------------


#cmake_minimum_required(VERSION 3.7)
#project(test_kalman_filters LANGUAGES CXX)
#message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}")

#enable_testing(true)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

#find_package (Boost COMPONENTS unit_test_framework REQUIRED)
#find_package(BLAS REQUIRED)
#find_package(Armadillo REQUIRED)
#include_directories(/usr/include/python3.8)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

#add_executable(${PROJECT_NAME}
#    test_kalman_filters.cpp
#    )
#add_test(NAME test_kalman_filters COMMAND test_kalman_filters)

#target_link_libraries(${PROJECT_NAME}
#    kalman_filters
#    libpython3.8.so
#    ${Boost_LIBRARIES}
#    )

#target_include_directories(${PROJECT_NAME}
#    PUBLIC
#        ${CMAKE_CURRENT_SOURCE_DIR}/../src/3d_party
#        ${CMAKE_CURRENT_SOURCE_DIR}/../src/compare
#    )
