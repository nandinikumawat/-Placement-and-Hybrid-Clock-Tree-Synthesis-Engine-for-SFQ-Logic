cmake_minimum_required(VERSION 3.24)

project(vda2_project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(FASTPLACE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/fastplace)
set(SFQPLACE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/sfqplace)

find_package(CGAL REQUIRED)

# Import hmetis static library 
if (NOT HMETIS_PATH)
    message(WARNING "HMETIS_PATH not definied (path to libhmetis.a). Will not link to HMETIS directly")
elseif(LINK_HMETIS)
    add_library(hmetis STATIC IMPORTED)
    set_target_properties(hmetis PROPERTIES IMPORTED_LOCATION ${HMETIS_PATH})
endif()

add_library(fastplace
    ${FASTPLACE_ROOT}/src/suraj_parser.cpp
    ${FASTPLACE_ROOT}/src/matrix.cpp
    ${FASTPLACE_ROOT}/src/placer.cpp
)
target_include_directories(fastplace PRIVATE ${FASTPLACE_ROOT}/include)

add_executable(PA3 ${FASTPLACE_ROOT}/src/main.cpp)
target_include_directories(PA3 PRIVATE ${FASTPLACE_ROOT}/include)
target_link_libraries(PA3 fastplace)

add_executable(sfqplace
    ${SFQPLACE_ROOT}/src/main.cpp
    ${SFQPLACE_ROOT}/src/grouping.cpp
    ${SFQPLACE_ROOT}/src/netlist.cpp
    ${SFQPLACE_ROOT}/src/partitioning.cpp
    ${SFQPLACE_ROOT}/src/supercells.cpp
)
target_include_directories(sfqplace PRIVATE
    ${SFQPLACE_ROOT}/include
    ${FASTPLACE_ROOT}/include
)
target_link_libraries(sfqplace fastplace gmp)
if (LINK_HMETIS)
    target_link_libraries(sfqplace hmetis)
    target_compile_definitions(sfqplace LINK_TO_HMETIS)
endif()

