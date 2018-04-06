cmake_minimum_required (VERSION 2.6)

# uncomment if on macOS with gcc installed
#set(CMAKE_C_COMPILER /usr/local/bin/gcc-7)
#set(CMAKE_CXX_COMPILER /usr/local/bin/g++-7)

project (Sig18)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pthread")

include_directories(include)

# add source files here!
add_library(database operators.cc relation.cc joiner.cc parser.cc worker.cc scheduler.cc)
target_include_directories(database PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

#OPTION(FORCE_TESTS "Build tests, regardless of build type." ON)
#if (CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]" OR FORCE_TESTS)
    #add_subdirectory(test)
#endif()

add_executable(Driver main.cc)
target_link_libraries(Driver database)

# Interactive command line tool to translate our query format to SQL queries
#add_executable(Query2SQL Query2SQL.cpp)
#target_link_libraries(Query2SQL database)

# Test harness
add_executable(harness harness.cc)

ADD_CUSTOM_TARGET(link_target ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/workloads
  ${CMAKE_CURRENT_BINARY_DIR}/workloads)
