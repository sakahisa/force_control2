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
include CMakeFiles/FastGetSetDCM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FastGetSetDCM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FastGetSetDCM.dir/flags.make

CMakeFiles/FastGetSetDCM.dir/fast.cpp.o: CMakeFiles/FastGetSetDCM.dir/flags.make
CMakeFiles/FastGetSetDCM.dir/fast.cpp.o: ../fast.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/FastGetSetDCM.dir/fast.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -o CMakeFiles/FastGetSetDCM.dir/fast.cpp.o -c /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/fast.cpp

CMakeFiles/FastGetSetDCM.dir/fast.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FastGetSetDCM.dir/fast.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -E /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/fast.cpp > CMakeFiles/FastGetSetDCM.dir/fast.cpp.i

CMakeFiles/FastGetSetDCM.dir/fast.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FastGetSetDCM.dir/fast.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -S /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/fast.cpp -o CMakeFiles/FastGetSetDCM.dir/fast.cpp.s

CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.requires:
.PHONY : CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.requires

CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.provides: CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.requires
	$(MAKE) -f CMakeFiles/FastGetSetDCM.dir/build.make CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.provides.build
.PHONY : CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.provides

CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.provides.build: CMakeFiles/FastGetSetDCM.dir/fast.cpp.o

CMakeFiles/FastGetSetDCM.dir/main.cpp.o: CMakeFiles/FastGetSetDCM.dir/flags.make
CMakeFiles/FastGetSetDCM.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/FastGetSetDCM.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -o CMakeFiles/FastGetSetDCM.dir/main.cpp.o -c /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/main.cpp

CMakeFiles/FastGetSetDCM.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FastGetSetDCM.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -E /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/main.cpp > CMakeFiles/FastGetSetDCM.dir/main.cpp.i

CMakeFiles/FastGetSetDCM.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FastGetSetDCM.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -S /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/main.cpp -o CMakeFiles/FastGetSetDCM.dir/main.cpp.s

CMakeFiles/FastGetSetDCM.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/FastGetSetDCM.dir/main.cpp.o.requires

CMakeFiles/FastGetSetDCM.dir/main.cpp.o.provides: CMakeFiles/FastGetSetDCM.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/FastGetSetDCM.dir/build.make CMakeFiles/FastGetSetDCM.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/FastGetSetDCM.dir/main.cpp.o.provides

CMakeFiles/FastGetSetDCM.dir/main.cpp.o.provides.build: CMakeFiles/FastGetSetDCM.dir/main.cpp.o

# Object files for target FastGetSetDCM
FastGetSetDCM_OBJECTS = \
"CMakeFiles/FastGetSetDCM.dir/fast.cpp.o" \
"CMakeFiles/FastGetSetDCM.dir/main.cpp.o"

# External object files for target FastGetSetDCM
FastGetSetDCM_EXTERNAL_OBJECTS =

sdk/lib/naoqi/libFastGetSetDCM.so: CMakeFiles/FastGetSetDCM.dir/fast.cpp.o
sdk/lib/naoqi/libFastGetSetDCM.so: CMakeFiles/FastGetSetDCM.dir/main.cpp.o
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalproxies.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalcommon.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalsoap.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/librttools.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalthread.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_signals-mt.a
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_program_options-mt.a
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalvalue.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libtinyxml.so
sdk/lib/naoqi/libFastGetSetDCM.so: /usr/lib/i386-linux-gnu/librt.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libqi.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_filesystem-mt.a
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_thread-mt.a
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libboost_system-mt.a
sdk/lib/naoqi/libFastGetSetDCM.so: /usr/lib/i386-linux-gnu/libdl.so
sdk/lib/naoqi/libFastGetSetDCM.so: /home/saka/Nao/naoqi-sdk-1.14.5-linux32/lib/libalerror.so
sdk/lib/naoqi/libFastGetSetDCM.so: CMakeFiles/FastGetSetDCM.dir/build.make
sdk/lib/naoqi/libFastGetSetDCM.so: CMakeFiles/FastGetSetDCM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library sdk/lib/naoqi/libFastGetSetDCM.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FastGetSetDCM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FastGetSetDCM.dir/build: sdk/lib/naoqi/libFastGetSetDCM.so
.PHONY : CMakeFiles/FastGetSetDCM.dir/build

CMakeFiles/FastGetSetDCM.dir/requires: CMakeFiles/FastGetSetDCM.dir/fast.cpp.o.requires
CMakeFiles/FastGetSetDCM.dir/requires: CMakeFiles/FastGetSetDCM.dir/main.cpp.o.requires
.PHONY : CMakeFiles/FastGetSetDCM.dir/requires

CMakeFiles/FastGetSetDCM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FastGetSetDCM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FastGetSetDCM.dir/clean

CMakeFiles/FastGetSetDCM.dir/depend:
	cd /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain /home/saka/Nao/naoqi-sdk-1.14.5-linux32/workspace/fast/build-crosstoolchain/CMakeFiles/FastGetSetDCM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FastGetSetDCM.dir/depend

