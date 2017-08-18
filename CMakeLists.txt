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
    MESSAGE(FATAL_ERROR "GCTA can't build on Windows with CMake currently, you can build with Visual Studio in gcta_win64 directory")
endif()

SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
IF(NOT EIGEN3_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "Specify environment variable EIGEN3_INCLUDE_DIR to the location of Eigen3 package")
ENDIF()
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")

SET(MKLROOT "$ENV{MKLROOT}")
IF(NOT MKLROOT)
    MESSAGE(FATAL_ERROR "Specify environment variable MKLROOT to the location of Intel MKL package")
ENDIF()
INCLUDE_DIRECTORIES("${MKLROOT}/include")

if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

find_library(lib_lp64 NAMES mkl_intel_lp64 PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
if(APPLE)
    find_library(lib_thread NAMES mkl_intel_thread PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
else()
    find_library(lib_thread NAMES mkl_gnu_thread PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
endif()
find_library(lib_core NAMES mkl_core PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
set(MKL_LIB ${lib_lp64} ${lib_thread} ${lib_core})

set(CMAKE_BUILD_TYPE RELEASE)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -fopenmp -msse4.2 -DMKL_LP64")
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -O3 -flto -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0 -g")
    # message(FALTAL_ERROR "GCTA has an known issue compile with Clang, we are checking the code")
    message(WARNING "GCTA has an known issue compile with Clang, we are checking the code")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -fopenmp -msse4.2 -pthread -DMKL_LP64 -fvisibility=hidden -fvisibility-inlines-hidden")
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -O3 -flto -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0 -Og")
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -axAVX2 -pthread")
    set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -O3 -ipo -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0 -Og")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(FALTAL_ERROR "GCTA can not compile by Microsft C++ compiler in CMake currently, please use gcta_win64 folder to build")
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

add_subdirectory(main)

if(APPLE)
    link_directories(${MKLROOT}/lib)
endif()

add_executable(gcta64 ${MAIN_SOURCE})

if(APPLE)
    target_link_libraries(gcta64 ${libs_list} mainV1 z ${MKL_LIB} iomp5 pthread m dl )
else()
    target_link_libraries(gcta64 ${libs_list} mainV1 -static z -Wl,--start-group ${MKL_LIB} -Wl,--end-group gomp pthread m dl )
endif()
#target_link_libraries(gcta64 ${libs_list} mainV1)

# Testing has some problems currently
#enable_testing()
#ADD_SUBDIRECTORY(test)

IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

