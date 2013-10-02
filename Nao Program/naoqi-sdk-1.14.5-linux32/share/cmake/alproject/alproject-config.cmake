# This is an autogenerated file. Do not edit

get_filename_component(_cur_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
set(_root_dir "${_cur_dir}/../../../")
get_filename_component(ROOT_DIR ${_root_dir} ABSOLUTE)

 
set(ALPROJECT_INCLUDE_DIRS "${ROOT_DIR}/include;" CACHE STRING "" FORCE)
mark_as_advanced(ALPROJECT_INCLUDE_DIRS)
   

find_library(ALPROJECT_DEBUG_LIBRARY alproject_d)
find_library(ALPROJECT_LIBRARY       alproject)


if (ALPROJECT_DEBUG_LIBRARY)
  set(ALPROJECT_LIBRARIES optimized;${ALPROJECT_LIBRARY};debug;${ALPROJECT_DEBUG_LIBRARY})
else()
  set(ALPROJECT_LIBRARIES ${ALPROJECT_LIBRARY})
endif()

set(ALPROJECT_LIBRARIES ${ALPROJECT_LIBRARIES} CACHE INTERNAL "" FORCE)
 
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ALPROJECT DEFAULT_MSG
  ALPROJECT_LIBRARIES
  ALPROJECT_INCLUDE_DIRS
)
set(ALPROJECT_PACKAGE_FOUND ${ALPROJECT_FOUND} CACHE INTERNAL "" FORCE)
 
set(ALPROJECT_DEPENDS "BOOST;ARCHIVE;ALFILE;QIPROJECT-QT;ALTOOLS;ALCORE;ALERROR;qiproject;qt_qtcore;alserial;qi;boost_signals;archive;QI;TINYXML;BOOST_FILESYSTEM;BOOST_THREAD;PTHREAD;DL;ZLIB;boost_filesystem;BOOST_SYSTEM" CACHE INTERNAL "" FORCE)
 