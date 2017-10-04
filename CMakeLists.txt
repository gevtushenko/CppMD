cmake_minimum_required(VERSION 3.2)
project(CppMD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

set(SOURCE_FILES
        includes/random/random.h
        sources/random/random.cpp
        includes/core/atom.h sources/core/atom.cpp sources/solvers/md_solver.cpp includes/solvers/md_solver.h)

include_directories(includes)

add_library(${CMAKE_PROJECT_NAME}_lib STATIC ${SOURCE_FILES})
add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES} main.cpp )
target_link_libraries(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}_lib)

find_package(Boost COMPONENTS unit_test_framework)

# Tests
if(Boost_FOUND)
    if(TEST_COV)
        SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage -std=c++11") # For code coverage
        SET(CMAKE_C_FLAGS "-g -O0 -Wall -W -fprofile-arcs -ftest-coverage -std=c++11")# For code coverage
        SET(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage -std=c++11")       # For code coverage
    endif()

    include_directories(${Boost_INCLUDE_DIR})
    file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tests/*.cpp)

    foreach(testSrc ${TEST_SRCS})
        get_filename_component(testName ${testSrc} NAME_WE)
        add_executable(${testName} ${testSrc})
        target_link_libraries(${testName} ${Boost_LIBRARIES} ${CMAKE_PROJECT_NAME}_lib)

        set_target_properties(
                ${testName} PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/tests/
        )
        add_test(
                NAME ${testName}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/tests/
                COMMAND ${PROJECT_SOURCE_DIR}/bin/tests/${testName}
        )
    endforeach()
    enable_testing()
endif()

# Docs
find_package(Doxygen)

if(${DOXYGEN_FOUND})
    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc)
endif()
