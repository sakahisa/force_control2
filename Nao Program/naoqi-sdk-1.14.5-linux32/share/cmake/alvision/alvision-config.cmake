# This is an autogenerated file. Do not edit

get_filename_component(_cur_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
set(_root_dir "${_cur_dir}/../../../")
get_filename_component(ROOT_DIR ${_root_dir} ABSOLUTE)

 
set(ALVISION_INCLUDE_DIRS "${ROOT_DIR}/include;" CACHE STRING "" FORCE)
mark_as_advanced(ALVISION_INCLUDE_DIRS)
   

find_library(ALVISION_DEBUG_LIBRARY alvision_d)
find_library(ALVISION_LIBRARY       alvision)


if (ALVISION_DEBUG_LIBRARY)
  set(ALVISION_LIBRARIES optimized;${ALVISION_LIBRARY};debug;${ALVISION_DEBUG_LIBRARY})
else()
  set(ALVISION_LIBRARIES ${ALVISION_LIBRARY})
endif()

set(ALVISION_LIBRARIES ${ALVISION_LIBRARIES} CACHE INTERNAL "" FORCE)
 
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ALVISION DEFAULT_MSG
  ALVISION_LIBRARIES
  ALVISION_INCLUDE_DIRS
)
set(ALVISION_PACKAGE_FOUND ${ALVISION_FOUND} CACHE INTERNAL "" FORCE)
 
set(ALVISION_DEPENDS "ALEXTRACTOR;ALPROXIES;ALCOMMON;ALCOMMON-INTERNAL;ALSOAP;RTTOOLS;ALTHREAD;BOOST_SIGNALS;BOOST_PROGRAM_OPTIONS;ALVALUE;TINYXML;RT;QI;BOOST_FILESYSTEM;BOOST_THREAD;BOOST_SYSTEM;PTHREAD;DL;ALERROR;BOOST" CACHE INTERNAL "" FORCE)
 