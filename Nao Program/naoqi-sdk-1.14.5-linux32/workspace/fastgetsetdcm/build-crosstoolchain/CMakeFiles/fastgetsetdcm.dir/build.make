# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain

# Include any dependencies generated for this target.
include CMakeFiles/fastgetsetdcm.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/fastgetsetdcm.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fastgetsetdcm.dir/flags.make

CMakeFiles/fastgetsetdcm.dir/main.cpp.o: CMakeFiles/fastgetsetdcm.dir/flags.make
CMakeFiles/fastgetsetdcm.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/fastgetsetdcm.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -o CMakeFiles/fastgetsetdcm.dir/main.cpp.o -c /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/main.cpp

CMakeFiles/fastgetsetdcm.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastgetsetdcm.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -E /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/main.cpp > CMakeFiles/fastgetsetdcm.dir/main.cpp.i

CMakeFiles/fastgetsetdcm.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastgetsetdcm.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -S /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/main.cpp -o CMakeFiles/fastgetsetdcm.dir/main.cpp.s

CMakeFiles/fastgetsetdcm.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/fastgetsetdcm.dir/main.cpp.o.requires

CMakeFiles/fastgetsetdcm.dir/main.cpp.o.provides: CMakeFiles/fastgetsetdcm.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/fastgetsetdcm.dir/build.make CMakeFiles/fastgetsetdcm.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/fastgetsetdcm.dir/main.cpp.o.provides

CMakeFiles/fastgetsetdcm.dir/main.cpp.o.provides.build: CMakeFiles/fastgetsetdcm.dir/main.cpp.o

CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o: CMakeFiles/fastgetsetdcm.dir/flags.make
CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o: ../fastgetsetdcm.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -o CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o -c /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/fastgetsetdcm.cpp

CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -E /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/fastgetsetdcm.cpp > CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.i

CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -S /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/fastgetsetdcm.cpp -o CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.s

CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.requires:
.PHONY : CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.requires

CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.provides: CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.requires
	$(MAKE) -f CMakeFiles/fastgetsetdcm.dir/build.make CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.provides.build
.PHONY : CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.provides

CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.provides.build: CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o

# Object files for target fastgetsetdcm
fastgetsetdcm_OBJECTS = \
"CMakeFiles/fastgetsetdcm.dir/main.cpp.o" \
"CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o"

# External object files for target fastgetsetdcm
fastgetsetdcm_EXTERNAL_OBJECTS =

sdk/lib/naoqi/libfastgetsetdcm.so: CMakeFiles/fastgetsetdcm.dir/main.cpp.o
sdk/lib/naoqi/libfastgetsetdcm.so: CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalmemoryfastaccess.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalproxies.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalproxies.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalcommon.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalsoap.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/librttools.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalthread.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_signals-mt.a
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_program_options-mt.a
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalvalue.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libtinyxml.so
sdk/lib/naoqi/libfastgetsetdcm.so: /usr/lib/i386-linux-gnu/librt.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libqi.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_filesystem-mt.a
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_thread-mt.a
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_system-mt.a
sdk/lib/naoqi/libfastgetsetdcm.so: /usr/lib/i386-linux-gnu/libdl.so
sdk/lib/naoqi/libfastgetsetdcm.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalerror.so
sdk/lib/naoqi/libfastgetsetdcm.so: CMakeFiles/fastgetsetdcm.dir/build.make
sdk/lib/naoqi/libfastgetsetdcm.so: CMakeFiles/fastgetsetdcm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library sdk/lib/naoqi/libfastgetsetdcm.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fastgetsetdcm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fastgetsetdcm.dir/build: sdk/lib/naoqi/libfastgetsetdcm.so
.PHONY : CMakeFiles/fastgetsetdcm.dir/build

CMakeFiles/fastgetsetdcm.dir/requires: CMakeFiles/fastgetsetdcm.dir/main.cpp.o.requires
CMakeFiles/fastgetsetdcm.dir/requires: CMakeFiles/fastgetsetdcm.dir/fastgetsetdcm.cpp.o.requires
.PHONY : CMakeFiles/fastgetsetdcm.dir/requires

CMakeFiles/fastgetsetdcm.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fastgetsetdcm.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fastgetsetdcm.dir/clean

CMakeFiles/fastgetsetdcm.dir/depend:
	cd /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fastgetset/build-crosstoolchain/CMakeFiles/fastgetsetdcm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fastgetsetdcm.dir/depend

