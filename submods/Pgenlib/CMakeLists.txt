SET(CMAKE_INCLUDE_CURRENT_DIR ON)
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
## don't add other PLINK GPL v3 codes, may have license issue
file(GLOB SOURCES_PLINK "${CMAKE_CURRENT_SOURCE_DIR}/../plink-ng/2.0/pgenlib*.cc")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../plink-ng/2.0/plink2_base.cc")
list(APPEND SOURCES "${SOURCES_PLINK}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../plink-ng/2.0")

add_library(Pgenlib STATIC ${SOURCES})
