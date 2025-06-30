###############################################################################
# Find zstd
#
# This sets the following variables:
# zstd_FOUND - True if zstd was found.
# zstd_INCLUDE_DIRS - Directories containing the zstd include files.
# zstd_LIBRARIES - Libraries needed to use zstd.

find_path(
    zstd_INCLUDE_DIR zstd.h
    PATHS
        /opt/local/include
        /usr/local/include
        /usr/include
    PATH_SUFFIXES zstd
)

find_library(
    zstd_LIBRARY
    NAMES zstd
    PATHS
        /opt/local/lib
        /user/local/lib
        /usr/lib
)

# Plural forms
set(zstd_INCLUDE_DIRS ${zstd_INCLUDE_DIR})
set(zstd_LIBRARIES ${zstd_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( zstd
  FOUND_VAR zstd_FOUND
  REQUIRED_VARS zstd_INCLUDE_DIR zstd_LIBRARY
)
