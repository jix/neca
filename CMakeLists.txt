cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(neca)
set(CMAKE_CXX_STANDARD 11)

find_package(GMP REQUIRED)
find_package(OpenMP)


add_executable(neca src/main.cpp)

target_include_directories(neca PRIVATE ${GMP_INCLUDE_DIRS})
target_link_libraries(neca ${GMP_LIBRARIES})

if(OpenMP_FOUND)
    target_compile_options(neca PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(neca ${OpenMP_CXX_LIBRARIES})
endif()
