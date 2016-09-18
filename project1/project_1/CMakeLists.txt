cmake_minimum_required(VERSION 2.8)

project(project_1)

find_package(Armadillo REQUIRED)
find_package(LAPACK    REQUIRED)
find_package(BLAS      REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11" )

# aux_source_directory(. SRC_LIST)
set(my_cpp_sources
    "${PROJECT_SOURCE_DIR}/main.cpp"
)

add_executable(${PROJECT_NAME} ${my_cpp_sources})
target_link_libraries(${PROJECT_NAME}
    ${ARMADILLO_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} )
