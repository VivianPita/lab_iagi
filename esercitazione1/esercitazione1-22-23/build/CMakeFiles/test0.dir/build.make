# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.29

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build

# Include any dependencies generated for this target.
include CMakeFiles/test0.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/test0.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/test0.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test0.dir/flags.make

CMakeFiles/test0.dir/src/test/test0.cpp.o: CMakeFiles/test0.dir/flags.make
CMakeFiles/test0.dir/src/test/test0.cpp.o: /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/src/test/test0.cpp
CMakeFiles/test0.dir/src/test/test0.cpp.o: CMakeFiles/test0.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test0.dir/src/test/test0.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test0.dir/src/test/test0.cpp.o -MF CMakeFiles/test0.dir/src/test/test0.cpp.o.d -o CMakeFiles/test0.dir/src/test/test0.cpp.o -c /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/src/test/test0.cpp

CMakeFiles/test0.dir/src/test/test0.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test0.dir/src/test/test0.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/src/test/test0.cpp > CMakeFiles/test0.dir/src/test/test0.cpp.i

CMakeFiles/test0.dir/src/test/test0.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test0.dir/src/test/test0.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/src/test/test0.cpp -o CMakeFiles/test0.dir/src/test/test0.cpp.s

# Object files for target test0
test0_OBJECTS = \
"CMakeFiles/test0.dir/src/test/test0.cpp.o"

# External object files for target test0
test0_EXTERNAL_OBJECTS =

/Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/test0: CMakeFiles/test0.dir/src/test/test0.cpp.o
/Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/test0: CMakeFiles/test0.dir/build.make
/Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/test0: libuwimg++.dylib
/Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/test0: CMakeFiles/test0.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/test0"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test0.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test0.dir/build: /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/test0
.PHONY : CMakeFiles/test0.dir/build

CMakeFiles/test0.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test0.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test0.dir/clean

CMakeFiles/test0.dir/depend:
	cd /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23 /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23 /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build /Users/alessandro/Desktop/esercitazione1/esercitazione1-22-23/build/CMakeFiles/test0.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/test0.dir/depend

