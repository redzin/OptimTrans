# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_SOURCE_DIR = "/home/jessica/johnnyredzin@gmail.com/MSc Thesis"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build"

# Include any dependencies generated for this target.
include CMakeFiles/main.o.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.o.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.o.dir/flags.make

CMakeFiles/main.o.dir/src/OTL.cpp.o: CMakeFiles/main.o.dir/flags.make
CMakeFiles/main.o.dir/src/OTL.cpp.o: ../src/OTL.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.o.dir/src/OTL.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.o.dir/src/OTL.cpp.o -c "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/src/OTL.cpp"

CMakeFiles/main.o.dir/src/OTL.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.o.dir/src/OTL.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/src/OTL.cpp" > CMakeFiles/main.o.dir/src/OTL.cpp.i

CMakeFiles/main.o.dir/src/OTL.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.o.dir/src/OTL.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/src/OTL.cpp" -o CMakeFiles/main.o.dir/src/OTL.cpp.s

CMakeFiles/main.o.dir/src/main.cpp.o: CMakeFiles/main.o.dir/flags.make
CMakeFiles/main.o.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.o.dir/src/main.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.o.dir/src/main.cpp.o -c "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/src/main.cpp"

CMakeFiles/main.o.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.o.dir/src/main.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/src/main.cpp" > CMakeFiles/main.o.dir/src/main.cpp.i

CMakeFiles/main.o.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.o.dir/src/main.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/src/main.cpp" -o CMakeFiles/main.o.dir/src/main.cpp.s

# Object files for target main.o
main_o_OBJECTS = \
"CMakeFiles/main.o.dir/src/OTL.cpp.o" \
"CMakeFiles/main.o.dir/src/main.cpp.o"

# External object files for target main.o
main_o_EXTERNAL_OBJECTS =

../main.o: CMakeFiles/main.o.dir/src/OTL.cpp.o
../main.o: CMakeFiles/main.o.dir/src/main.cpp.o
../main.o: CMakeFiles/main.o.dir/build.make
../main.o: /usr/lib/x86_64-linux-gnu/libtiff.so
../main.o: /usr/lib/x86_64-linux-gnu/libjpeg.so
../main.o: /usr/lib/x86_64-linux-gnu/libz.so
../main.o: /usr/lib/x86_64-linux-gnu/libpng.so
../main.o: /usr/lib/x86_64-linux-gnu/libz.so
../main.o: /usr/lib/x86_64-linux-gnu/libSM.so
../main.o: /usr/lib/x86_64-linux-gnu/libICE.so
../main.o: /usr/lib/x86_64-linux-gnu/libX11.so
../main.o: /usr/lib/x86_64-linux-gnu/libXext.so
../main.o: /usr/lib/x86_64-linux-gnu/liblapack.so
../main.o: /usr/lib/x86_64-linux-gnu/libblas.so
../main.o: /usr/lib/x86_64-linux-gnu/libblas.so
../main.o: /usr/lib/x86_64-linux-gnu/libpng.so
../main.o: /usr/lib/x86_64-linux-gnu/libSM.so
../main.o: /usr/lib/x86_64-linux-gnu/libICE.so
../main.o: /usr/lib/x86_64-linux-gnu/libX11.so
../main.o: /usr/lib/x86_64-linux-gnu/libXext.so
../main.o: /usr/lib/x86_64-linux-gnu/liblapack.so
../main.o: /usr/lib/x86_64-linux-gnu/libblas.so
../main.o: CMakeFiles/main.o.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../main.o"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.o.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.o.dir/build: ../main.o

.PHONY : CMakeFiles/main.o.dir/build

CMakeFiles/main.o.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.o.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.o.dir/clean

CMakeFiles/main.o.dir/depend:
	cd "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/jessica/johnnyredzin@gmail.com/MSc Thesis" "/home/jessica/johnnyredzin@gmail.com/MSc Thesis" "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build" "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build" "/home/jessica/johnnyredzin@gmail.com/MSc Thesis/build/CMakeFiles/main.o.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/main.o.dir/depend

