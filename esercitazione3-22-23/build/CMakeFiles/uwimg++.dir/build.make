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
CMAKE_SOURCE_DIR = /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build

# Include any dependencies generated for this target.
include CMakeFiles/uwimg++.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/uwimg++.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/uwimg++.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/uwimg++.dir/flags.make

CMakeFiles/uwimg++.dir/src/utils.cpp.o: CMakeFiles/uwimg++.dir/flags.make
CMakeFiles/uwimg++.dir/src/utils.cpp.o: /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/utils.cpp
CMakeFiles/uwimg++.dir/src/utils.cpp.o: CMakeFiles/uwimg++.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/uwimg++.dir/src/utils.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/uwimg++.dir/src/utils.cpp.o -MF CMakeFiles/uwimg++.dir/src/utils.cpp.o.d -o CMakeFiles/uwimg++.dir/src/utils.cpp.o -c /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/utils.cpp

CMakeFiles/uwimg++.dir/src/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/uwimg++.dir/src/utils.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/utils.cpp > CMakeFiles/uwimg++.dir/src/utils.cpp.i

CMakeFiles/uwimg++.dir/src/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/uwimg++.dir/src/utils.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/utils.cpp -o CMakeFiles/uwimg++.dir/src/utils.cpp.s

CMakeFiles/uwimg++.dir/src/load_image.cpp.o: CMakeFiles/uwimg++.dir/flags.make
CMakeFiles/uwimg++.dir/src/load_image.cpp.o: /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/load_image.cpp
CMakeFiles/uwimg++.dir/src/load_image.cpp.o: CMakeFiles/uwimg++.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/uwimg++.dir/src/load_image.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/uwimg++.dir/src/load_image.cpp.o -MF CMakeFiles/uwimg++.dir/src/load_image.cpp.o.d -o CMakeFiles/uwimg++.dir/src/load_image.cpp.o -c /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/load_image.cpp

CMakeFiles/uwimg++.dir/src/load_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/uwimg++.dir/src/load_image.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/load_image.cpp > CMakeFiles/uwimg++.dir/src/load_image.cpp.i

CMakeFiles/uwimg++.dir/src/load_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/uwimg++.dir/src/load_image.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/load_image.cpp -o CMakeFiles/uwimg++.dir/src/load_image.cpp.s

CMakeFiles/uwimg++.dir/src/process_image.cpp.o: CMakeFiles/uwimg++.dir/flags.make
CMakeFiles/uwimg++.dir/src/process_image.cpp.o: /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/process_image.cpp
CMakeFiles/uwimg++.dir/src/process_image.cpp.o: CMakeFiles/uwimg++.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/uwimg++.dir/src/process_image.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/uwimg++.dir/src/process_image.cpp.o -MF CMakeFiles/uwimg++.dir/src/process_image.cpp.o.d -o CMakeFiles/uwimg++.dir/src/process_image.cpp.o -c /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/process_image.cpp

CMakeFiles/uwimg++.dir/src/process_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/uwimg++.dir/src/process_image.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/process_image.cpp > CMakeFiles/uwimg++.dir/src/process_image.cpp.i

CMakeFiles/uwimg++.dir/src/process_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/uwimg++.dir/src/process_image.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/process_image.cpp -o CMakeFiles/uwimg++.dir/src/process_image.cpp.s

CMakeFiles/uwimg++.dir/src/filter_image.cpp.o: CMakeFiles/uwimg++.dir/flags.make
CMakeFiles/uwimg++.dir/src/filter_image.cpp.o: /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/filter_image.cpp
CMakeFiles/uwimg++.dir/src/filter_image.cpp.o: CMakeFiles/uwimg++.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/uwimg++.dir/src/filter_image.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/uwimg++.dir/src/filter_image.cpp.o -MF CMakeFiles/uwimg++.dir/src/filter_image.cpp.o.d -o CMakeFiles/uwimg++.dir/src/filter_image.cpp.o -c /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/filter_image.cpp

CMakeFiles/uwimg++.dir/src/filter_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/uwimg++.dir/src/filter_image.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/filter_image.cpp > CMakeFiles/uwimg++.dir/src/filter_image.cpp.i

CMakeFiles/uwimg++.dir/src/filter_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/uwimg++.dir/src/filter_image.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/filter_image.cpp -o CMakeFiles/uwimg++.dir/src/filter_image.cpp.s

CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o: CMakeFiles/uwimg++.dir/flags.make
CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o: /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/edge_detection.cpp
CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o: CMakeFiles/uwimg++.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o -MF CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o.d -o CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o -c /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/edge_detection.cpp

CMakeFiles/uwimg++.dir/src/edge_detection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/uwimg++.dir/src/edge_detection.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/edge_detection.cpp > CMakeFiles/uwimg++.dir/src/edge_detection.cpp.i

CMakeFiles/uwimg++.dir/src/edge_detection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/uwimg++.dir/src/edge_detection.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/src/edge_detection.cpp -o CMakeFiles/uwimg++.dir/src/edge_detection.cpp.s

# Object files for target uwimg++
uwimg_______OBJECTS = \
"CMakeFiles/uwimg++.dir/src/utils.cpp.o" \
"CMakeFiles/uwimg++.dir/src/load_image.cpp.o" \
"CMakeFiles/uwimg++.dir/src/process_image.cpp.o" \
"CMakeFiles/uwimg++.dir/src/filter_image.cpp.o" \
"CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o"

# External object files for target uwimg++
uwimg_______EXTERNAL_OBJECTS =

libuwimg++.dylib: CMakeFiles/uwimg++.dir/src/utils.cpp.o
libuwimg++.dylib: CMakeFiles/uwimg++.dir/src/load_image.cpp.o
libuwimg++.dylib: CMakeFiles/uwimg++.dir/src/process_image.cpp.o
libuwimg++.dylib: CMakeFiles/uwimg++.dir/src/filter_image.cpp.o
libuwimg++.dylib: CMakeFiles/uwimg++.dir/src/edge_detection.cpp.o
libuwimg++.dylib: CMakeFiles/uwimg++.dir/build.make
libuwimg++.dylib: CMakeFiles/uwimg++.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX shared library libuwimg++.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/uwimg++.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/uwimg++.dir/build: libuwimg++.dylib
.PHONY : CMakeFiles/uwimg++.dir/build

CMakeFiles/uwimg++.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/uwimg++.dir/cmake_clean.cmake
.PHONY : CMakeFiles/uwimg++.dir/clean

CMakeFiles/uwimg++.dir/depend:
	cd /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23 /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23 /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build /Users/alessandro/Desktop/esercitazione3/esercitazione3-22-23/build/CMakeFiles/uwimg++.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/uwimg++.dir/depend

