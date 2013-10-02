## Copyright (c) 2012 Aldebaran Robotics. All rights reserved.
## Use of this source code is governed by a BSD-style license that can be
## found in the COPYING file.

get_filename_component(_ROOT_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
include("${_ROOT_DIR}/qtutils.cmake")

set(_suffix "QTMULTIMEDIA")
set(_libame "QtMultimedia")

qt_flib(${_suffix} ${_libame})
if(UNIX AND NOT APPLE)
  qi_set_global(QT_QTMULTIMEDIA_DEPENDS "ALSALIB")
endif()
export_lib(QT_${_suffix})
set(_ROOT_DIR)
