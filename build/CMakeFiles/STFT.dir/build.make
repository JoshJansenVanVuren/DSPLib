# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/josh/devel/DSPLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/josh/devel/DSPLib/build

# Include any dependencies generated for this target.
include CMakeFiles/STFT.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/STFT.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/STFT.dir/flags.make

CMakeFiles/STFT.dir/lib/stft.cpp.o: CMakeFiles/STFT.dir/flags.make
CMakeFiles/STFT.dir/lib/stft.cpp.o: ../lib/stft.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/devel/DSPLib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/STFT.dir/lib/stft.cpp.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/STFT.dir/lib/stft.cpp.o -c /home/josh/devel/DSPLib/lib/stft.cpp

CMakeFiles/STFT.dir/lib/stft.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/STFT.dir/lib/stft.cpp.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/devel/DSPLib/lib/stft.cpp > CMakeFiles/STFT.dir/lib/stft.cpp.i

CMakeFiles/STFT.dir/lib/stft.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/STFT.dir/lib/stft.cpp.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/devel/DSPLib/lib/stft.cpp -o CMakeFiles/STFT.dir/lib/stft.cpp.s

CMakeFiles/STFT.dir/lib/stft.cpp.o.requires:

.PHONY : CMakeFiles/STFT.dir/lib/stft.cpp.o.requires

CMakeFiles/STFT.dir/lib/stft.cpp.o.provides: CMakeFiles/STFT.dir/lib/stft.cpp.o.requires
	$(MAKE) -f CMakeFiles/STFT.dir/build.make CMakeFiles/STFT.dir/lib/stft.cpp.o.provides.build
.PHONY : CMakeFiles/STFT.dir/lib/stft.cpp.o.provides

CMakeFiles/STFT.dir/lib/stft.cpp.o.provides.build: CMakeFiles/STFT.dir/lib/stft.cpp.o


# Object files for target STFT
STFT_OBJECTS = \
"CMakeFiles/STFT.dir/lib/stft.cpp.o"

# External object files for target STFT
STFT_EXTERNAL_OBJECTS =

../lib/libSTFT.a: CMakeFiles/STFT.dir/lib/stft.cpp.o
../lib/libSTFT.a: CMakeFiles/STFT.dir/build.make
../lib/libSTFT.a: CMakeFiles/STFT.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/devel/DSPLib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ../lib/libSTFT.a"
	$(CMAKE_COMMAND) -P CMakeFiles/STFT.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/STFT.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/STFT.dir/build: ../lib/libSTFT.a

.PHONY : CMakeFiles/STFT.dir/build

CMakeFiles/STFT.dir/requires: CMakeFiles/STFT.dir/lib/stft.cpp.o.requires

.PHONY : CMakeFiles/STFT.dir/requires

CMakeFiles/STFT.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/STFT.dir/cmake_clean.cmake
.PHONY : CMakeFiles/STFT.dir/clean

CMakeFiles/STFT.dir/depend:
	cd /home/josh/devel/DSPLib/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/devel/DSPLib /home/josh/devel/DSPLib /home/josh/devel/DSPLib/build /home/josh/devel/DSPLib/build /home/josh/devel/DSPLib/build/CMakeFiles/STFT.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/STFT.dir/depend

