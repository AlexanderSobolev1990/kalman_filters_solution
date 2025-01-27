cmake_minimum_required(VERSION 3.7)
project(kalman_filters_solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE) # Если не указан тип сборки, собрать в режиме Release
  set(CMAKE_BUILD_TYPE Release)
endif()
# message(WARNING "Build mode: ${CMAKE_BUILD_TYPE}")
message("Build mode: ${CMAKE_BUILD_TYPE}")

# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/lib)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease") # Актуально для Linux
    option(BUILD_DOCUMENTATION "Build project documentation (Requires Doxygen)" ON) # Опционально построение документации: ON|OFF
    # SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp -O3")
    # SET(GCC_COVERAGE_LINK_FLAGS "-fopenmp")    
    SET(GCC_COVERAGE_COMPILE_FLAGS "-O3")       
    SET(CPACK_STRIP_FILES ON)
endif()

# add_compile_options(
#     # -Wall    
#     "$<$<CONFIG:RELEASE>:-O3>"
#     "$<$<CONFIG:DEBUG>:-O0>"
# )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

# Директивы
add_definitions(    
    -DNUM_OF_THREADS=4 # Число потоков распараллеливания вычислений
)

#-------------------------------------------------------------------------------
# Проекты
add_subdirectory(src/kalman_filters)
add_subdirectory(src/compare)

# Генерация документации
if(BUILD_DOCUMENTATION)
    find_package(Doxygen
        REQUIRED dot) # dot - требуется Graphviz dot утилита для построения графиков в документации
    if(DOXYGEN_FOUND)
        add_subdirectory(doc) # В поддиректории проекта /doc должен быть CMakeLists.txt с настройками Doxygen
    else()
        message(WARNING "Doxygen need to be installed to generate the doxygen documentation")
    endif()
endif()
