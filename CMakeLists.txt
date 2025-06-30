cmake_minimum_required(VERSION 3.16)
project(gcta)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# 设置调试模式
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif() 

if(WIN32 AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "Only the Clang compiler is supported on Windows, Please use Clang to compile this project.")
endif()



if (APPLE)
    # todo: link openblas for macOS
else()
    if(UNIX)
        set(MKL_THREADING "gnu_thread")
    endif()
    set(MKL_LINK "static")
    set(MKL_INTERFACE "lp64")
    find_package(MKL CONFIG REQUIRED QUIET) 
    if(MKL_FOUND)
        message(STATUS "MKL_INCLUDE: ${MKL_INCLUDE}")
        include_directories(${MKL_INCLUDE})
    endif()
endif()

find_package(Boost CONFIG REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(Eigen3 CONFIG REQUIRED)
if(Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIRS})
    message(STATUS "EIGEN3_INCLUDE_DIRS: ${EIGEN3_INCLUDE_DIRS}")
endif()

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    message(STATUS "ZLIB_INCLUDE_DIRS: ${ZLIB_INCLUDE_DIRS}")
    message(STATUS "ZLIB_LIBRARIES: ${ZLIB_LIBRARIES}")
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

find_package(SQLite3 REQUIRED)
if(SQLite3_FOUND)
    message(STATUS "SQLite3_INCLUDE_DIRS: ${SQLite3_INCLUDE_DIRS}")
    message(STATUS "SQLite3_LIBRARIES: ${SQLite3_LIBRARIES}")
    include_directories(${SQLite3_INCLUDE_DIRS})
endif()

find_package(GSL REQUIRED)
if(GSL_FOUND)
    message(STATUS "GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIRS}")
    message(STATUS "GSL_LIBRARIES: ${GSL_LIBRARIES}")
    include_directories(${GSL_INCLUDE_DIRS})
endif()

find_package(zstd REQUIRED)
if(zstd_FOUND)
    message(STATUS "zstd_INCLUDE_DIRS: ${zstd_INCLUDE_DIRS}")
    message(STATUS "zstd_LIBRARIES: ${zstd_LIBRARIES}")
    include_directories(${zstd_INCLUDE_DIRS})
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src/Pgenlib
    ${PROJECT_SOURCE_DIR}/thirdparty/plink-ng/2.0/include
)

add_subdirectory(src)

add_executable(${PROJECT_NAME} "")

set(MAIN_SOURCE "${PROJECT_SOURCE_DIR}/src/core/main.cpp")
file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/core/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")
list(APPEND SRCS "${PROJECT_SOURCE_DIR}/src/core/acat.c")

set(libs_list "")
foreach(_lib ${SRCS})
    get_filename_component(file_name "${_lib}" NAME_WE)
    string(TOLOWER ${file_name} lib_name)
    add_library(${lib_name} STATIC "${_lib}")
    list(APPEND libs_list ${lib_name})
endforeach(_lib)


target_link_libraries(${PROJECT_NAME} 
    mainV1
    ${libs_list} 
    Pgenlib
    ${ZLIB_LIBRARIES}
    ${SQLite3_LIBRARIES}
    ${zstd_LIBRARIES}
    ${GSL_LIBRARIES}
    MKL::MKL
)

target_sources(${PROJECT_NAME} 
  PRIVATE
    ${MAIN_SOURCE}
)