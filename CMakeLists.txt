CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)
SET(PROJECT_NAME_STR GCTA2)
PROJECT(${PROJECT_NAME_STR} C CXX)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_STANDARD 11)

# get system architecture
EXECUTE_PROCESS(COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE)
set(CMAKE_LIBRARY_ARCHITECTURE ${ARCHITECTURE})

# fix win stack size issue
if(WIN32)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--stack,10000000")
    MESSAGE(FATAL_ERROR "Use the Visual Studio project file which located in gcta_win64 directory")
endif()

if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux" OR "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    MESSAGE("Linux and Mac are supported.")
else()
    MESSAGE("Your system: ${CMAKE_SYSTEM_NAME}")
    MESSAGE(FATAL_ERROR "Only Linux and Mac are currently supported.")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    MESSAGE("Compiling by GCC${CMAKE_CXX_COMPILER_VERSION}")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0))
        MESSAGE(FATAL_ERROR "GCTA can't build without GCC version >= 6.1")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("Compiling by Clang${CMAKE_CXX_COMPILER_VERSION}")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0))
        MESSAGE(FATAL_ERROR "GCTA can't build without Clang version >= 8.0")
    endif()
else()
    MESSAGE(FATAL_ERROR "GCTA can only compile with GCC or Clang")
endif()

SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
IF(NOT EIGEN3_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "Specify environment variable EIGEN3_INCLUDE_DIR to the location of Eigen3 package")
ENDIF()
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")

SET(SPECTRA_LIB "$ENV{SPECTRA_LIB}")
IF(NOT SPECTRA_LIB)
    MESSAGE(FATAL_ERROR "Specify environment variable SPECTRA_LIB to the location of Spectra package")
ENDIF()
INCLUDE_DIRECTORIES("${SPECTRA_LIB}")


SET(BOOST_LIB "$ENV{BOOST_LIB}")
IF(NOT BOOST_LIB)
    MESSAGE(FATAL_ERROR "Specify environment variable BOOST_LIB to the location of BOOST library")
ENDIF()
INCLUDE_DIRECTORIES("${BOOST_LIB}")


# PLINK2.0 pgen library;
INCLUDE_DIRECTORIES("submods/plink-ng/2.0")
INCLUDE_DIRECTORIES("submods/plink-ng/2.0/simde")


message(STATUS "ARCHITECTURE:${CMAKE_LIBRARY_ARCHITECTURE}")
if(CMAKE_LIBRARY_ARCHITECTURE STREQUAL "aarch64")
    SET(OPENBLAS "$ENV{OPENBLAS}")
    IF(NOT OPENBLAS)
        MESSAGE(FATAL_ERROR "Specify environment variable OPENBLAS to the location of OPENBLAS package")
    ENDIF()
    INCLUDE_DIRECTORIES("${OPENBLAS}/include")

    find_library(lib_openblas NAMES openblas PATHS "${OPENBLAS}" "${OPENBLAS}/lib")
    set(BLAS_LIB ${lib_openblas})

    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
        add_compile_options(-fopenmp)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_implicit_dylibs")
        if(NOT APPLE)
            # don't delete this it will result in _implibs
            set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -pthread -fvisibility=hidden -fvisibility-inlines-hidden")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG") #-flto
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fno-inline -fno-implicit-inline-templates -g3")
        add_compile_options(-fopenmp)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ipo -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
    endif()

else()
    SET(MKLROOT "$ENV{MKLROOT}")
    IF(NOT MKLROOT)
        MESSAGE(FATAL_ERROR "Specify environment variable MKLROOT to the location of Intel MKL package")
    ENDIF()
    INCLUDE_DIRECTORIES("${MKLROOT}/include")

    find_library(lib_lp64 NAMES mkl_intel_lp64 PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR APPLE)
        find_library(lib_thread NAMES mkl_intel_thread PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
    else()
        find_library(lib_thread NAMES mkl_gnu_thread PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
    endif()
    find_library(lib_core NAMES mkl_core PATHS "${MKLROOT}/lib/intel64" "${MKLROOT}/lib")
    set(BLAS_LIB ${lib_lp64} ${lib_thread} ${lib_core})

    IF(NOT SSE)
        set(SSE "-msse2")
        message("SSE build")
    ENDIF()

    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ${SSE} -DMKL_LP64 -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
        add_compile_options(-fopenmp)
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_implicit_dylibs")
        if(NOT APPLE)
            # don't delete this it will result in _implibs
            set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 ${SSE} -pthread -DMKL_LP64 -fvisibility=hidden -fvisibility-inlines-hidden")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG") #-flto
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fno-inline -fno-implicit-inline-templates -g3")
        add_compile_options(-fopenmp)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -axAVX2 -pthread")
        set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ipo -DNDEBUG")
        set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        message(FALTAL_ERROR "GCTA can not compile by Microsft C++ compiler in CMake currently, please use gcta_win64 folder to build")
    endif()
endif()

IF(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
    message("Setting to release build")
ENDIF()

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
set(MAIN_SOURCE "${PROJECT_SOURCE_DIR}/src/main.cpp")
include_directories(${COMMON_INCLUDES})


file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")
list(APPEND SRCS "${PROJECT_SOURCE_DIR}/src/acat.c")

set(libs_list "")
foreach(_lib ${SRCS})
    get_filename_component(file_name "${_lib}" NAME_WE)
    string(TOLOWER ${file_name} lib_name)
    add_library(${lib_name} STATIC "${_lib}")
    list(APPEND libs_list ${lib_name})
endforeach(_lib)

add_subdirectory(main)
add_subdirectory(submods)

if(APPLE)
    link_directories(${MKLROOT}/lib)
    link_directories(/usr/local/lib)
else(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    #linux
    link_directories(${MKLROOT}/../../../lib/intel64)
endif()

add_executable(gcta64 ${MAIN_SOURCE})

if(APPLE)
    # iomp5
    target_link_libraries(gcta64 ${libs_list} Pgenlib mainV1 z ${BLAS_LIB} sqlite3 zstd_s omp pthread m dl )
else()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # error keeps only
        #target_link_libraries(gcta64 ${libs_list} mainV1 z ${MKL_LIB} iomp5 pthread m dl )
        #target_link_libraries(gcta64 mainV1 ${libs_list} z ${MKL_LIB} pthread m dl)
        # Dynamic link works
        #target_link_libraries(gcta64 mainV1 ${libs_list} -Wl,-Bstatic z -Wl,--start-group ${MKL_LIB} -Wl,--end-group iomp5 -Wl,-Bdynamic pthread m dl)
        target_link_libraries(gcta64 mainV1 ${libs_list} Pgenlib z sqlite3 zstd -Wl,--start-group ${BLAS_LIB} -Wl,--end-group iomp5 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive m dl)
    else()
        #target_link_libraries(gcta64 mainV1 ${libs_list} -Wl,-Bstatic z -Wl,--start-group ${MKL_LIB} -Wl,--end-group gomp -Wl,-Bdynamic pthread m dl)
        #target_link_libraries(gcta64 mainV1 ${libs_list} Pgenlib -static z sqlite3 zstd -Wl,--start-group ${MKL_LIB} -Wl,--end-group gomp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive m dl)
        target_link_libraries(gcta64 mainV1 ${libs_list} Pgenlib z sqlite3 zstd gsl gslcblas m -Wl,--start-group ${BLAS_LIB} -Wl,--end-group gomp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive m dl)
    endif()
endif()

# Testing has some problems currently
#enable_testing()
#ADD_SUBDIRECTORY(test)

IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

