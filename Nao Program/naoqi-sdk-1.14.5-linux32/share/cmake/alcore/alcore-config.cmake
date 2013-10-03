# This is an autogenerated file. Do not edit

get_filename_component(_cur_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
set(_root_dir "${_cur_dir}/../../../")
get_filename_component(ROOT_DIR ${_root_dir} ABSOLUTE)

 
set(ALCORE_INCLUDE_DIRS "${ROOT_DIR}/include;" CACHE STRING "" FORCE)
mark_as_advanced(ALCORE_INCLUDE_DIRS)
   
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ALCORE DEFAULT_MSG
  ALCORE_INCLUDE_DIRS
)
set(ALCORE_PACKAGE_FOUND ${ALCORE_FOUND} CACHE INTERNAL "" FORCE)
 
set(ALCORE_DEPENDS "QI;ALERROR" CACHE INTERNAL "" FORCE)
 
if(NOT ALCORE_I_KNOW_IT_IS_DEPRECATED AND QI_WARN_DEPRECATED)
  message(STATUS "
    ${CMAKE_CURRENT_SOURCE_DIR}: ALCORE is deprecated
    please use libqi.
  "
  )
  set(ALCORE_I_KNOW_IT_IS_DEPRECATED ON)
endif()