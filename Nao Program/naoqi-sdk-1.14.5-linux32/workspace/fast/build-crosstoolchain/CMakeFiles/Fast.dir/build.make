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
CMAKE_SOURCE_DIR = /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain

# Include any dependencies generated for this target.
include CMakeFiles/Fast.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Fast.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Fast.dir/flags.make

CMakeFiles/Fast.dir/fast.cpp.o: CMakeFiles/Fast.dir/flags.make
CMakeFiles/Fast.dir/fast.cpp.o: ../fast.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Fast.dir/fast.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -o CMakeFiles/Fast.dir/fast.cpp.o -c /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/fast.cpp

CMakeFiles/Fast.dir/fast.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Fast.dir/fast.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -E /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/fast.cpp > CMakeFiles/Fast.dir/fast.cpp.i

CMakeFiles/Fast.dir/fast.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Fast.dir/fast.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -S /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/fast.cpp -o CMakeFiles/Fast.dir/fast.cpp.s

CMakeFiles/Fast.dir/fast.cpp.o.requires:
.PHONY : CMakeFiles/Fast.dir/fast.cpp.o.requires

CMakeFiles/Fast.dir/fast.cpp.o.provides: CMakeFiles/Fast.dir/fast.cpp.o.requires
	$(MAKE) -f CMakeFiles/Fast.dir/build.make CMakeFiles/Fast.dir/fast.cpp.o.provides.build
.PHONY : CMakeFiles/Fast.dir/fast.cpp.o.provides

CMakeFiles/Fast.dir/fast.cpp.o.provides.build: CMakeFiles/Fast.dir/fast.cpp.o

CMakeFiles/Fast.dir/main.cpp.o: CMakeFiles/Fast.dir/flags.make
CMakeFiles/Fast.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Fast.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -o CMakeFiles/Fast.dir/main.cpp.o -c /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/main.cpp

CMakeFiles/Fast.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Fast.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -E /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/main.cpp > CMakeFiles/Fast.dir/main.cpp.i

CMakeFiles/Fast.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Fast.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -S /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/main.cpp -o CMakeFiles/Fast.dir/main.cpp.s

CMakeFiles/Fast.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/Fast.dir/main.cpp.o.requires

CMakeFiles/Fast.dir/main.cpp.o.provides: CMakeFiles/Fast.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/Fast.dir/build.make CMakeFiles/Fast.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/Fast.dir/main.cpp.o.provides

CMakeFiles/Fast.dir/main.cpp.o.provides.build: CMakeFiles/Fast.dir/main.cpp.o

# Object files for target Fast
Fast_OBJECTS = \
"CMakeFiles/Fast.dir/fast.cpp.o" \
"CMakeFiles/Fast.dir/main.cpp.o"

# External object files for target Fast
Fast_EXTERNAL_OBJECTS =

sdk/lib/naoqi/libFast.so: CMakeFiles/Fast.dir/fast.cpp.o
sdk/lib/naoqi/libFast.so: CMakeFiles/Fast.dir/main.cpp.o
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalproxies.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalcommon.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalsoap.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/librttools.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalthread.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_signals-mt.a
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_program_options-mt.a
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalvalue.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libtinyxml.so
sdk/lib/naoqi/libFast.so: /usr/lib/i386-linux-gnu/librt.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libqi.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_filesystem-mt.a
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_thread-mt.a
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_system-mt.a
sdk/lib/naoqi/libFast.so: /usr/lib/i386-linux-gnu/libdl.so
sdk/lib/naoqi/libFast.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalerror.so
sdk/lib/naoqi/libFast.so: CMakeFiles/Fast.dir/build.make
sdk/lib/naoqi/libFast.so: CMakeFiles/Fast.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library sdk/lib/naoqi/libFast.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Fast.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Fast.dir/build: sdk/lib/naoqi/libFast.so
.PHONY : CMakeFiles/Fast.dir/build

CMakeFiles/Fast.dir/requires: CMakeFiles/Fast.dir/fast.cpp.o.requires
CMakeFiles/Fast.dir/requires: CMakeFiles/Fast.dir/main.cpp.o.requires
.PHONY : CMakeFiles/Fast.dir/requires

CMakeFiles/Fast.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Fast.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Fast.dir/clean

CMakeFiles/Fast.dir/depend:
	cd /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain/CMakeFiles/Fast.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Fast.dir/depend
