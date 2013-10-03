# This is an autogenerated file. Do not edit

get_filename_component(_cur_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
set(_root_dir "${_cur_dir}/../../../")
get_filename_component(ROOT_DIR ${_root_dir} ABSOLUTE)

 
set(LIBLAUNCHER_INCLUDE_DIRS "${ROOT_DIR}/include;" CACHE STRING "" FORCE)
mark_as_advanced(ALLAUNCHER_INCLUDE_DIRS)
   

find_library(LIBLAUNCHER_DEBUG_LIBRARY allauncher_d)
find_library(LIBLAUNCHER_LIBRARY       allauncher)


if (LIBLAUNCHER_DEBUG_LIBRARY)
  set(LIBLAUNCHER_LIBRARIES optimized;${LIBLAUNCHER_LIBRARY};debug;${LIBLAUNCHER_DEBUG_LIBRARY})
else()
  set(LIBLAUNCHER_LIBRARIES ${LIBLAUNCHER_LIBRARY})
endif()

set(LIBLAUNCHER_LIBRARIES ${LIBLAUNCHER_LIBRARIES} CACHE INTERNAL "" FORCE)
 
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBLAUNCHER DEFAULT_MSG
  LIBLAUNCHER_LIBRARIES
  LIBLAUNCHER_INCLUDE_DIRS
)
set(LIBLAUNCHER_PACKAGE_FOUND ${LIBLAUNCHER_FOUND} CACHE INTERNAL "" FORCE)
 
set(LIBLAUNCHER_DEPENDS "" CACHE INTERNAL "" FORCE)
 