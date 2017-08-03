CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)
SET(PROJECT_NAME_STR GCTA2)
PROJECT(${PROJECT_NAME_STR} C CXX)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_STANDARD 11)

# fix win stack size issue
if(WIN32)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--stack,10000000")
endif()

set(CMAKE_BUILD_TYPE RELEASE)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -msse4.2")
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -O3 -flto")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0 -g")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -msse4.2 -pthread")
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -O3 -flto")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0 -Og")
    #set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    #set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -axAVX2 -pthread")
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -O3 -ipo")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0 -Og")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(FALTAL_ERROR "GCTA can not compile by Microsft C++ compiler currently, we are debuging this")
endif()

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(MAIN_SOURCE "${PROJECT_SOURCE_DIR}/src/main.cpp")
include_directories(${COMMON_INCLUDES})

file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")

set(libs_list "")
foreach(_lib ${SRCS})
    get_filename_component(file_name "${_lib}" NAME_WE)
    string(TOLOWER ${file_name} lib_name)
    add_library(${lib_name} STATIC "${_lib}")
    list(APPEND libs_list ${lib_name})
endforeach(_lib)

add_executable(gcta2 ${MAIN_SOURCE})
target_link_libraries(gcta2 ${libs_list})

enable_testing()
ADD_SUBDIRECTORY(test)

IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

