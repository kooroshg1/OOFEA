# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /Users/gobalk/Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Users/gobalk/Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build

# Include any dependencies generated for this target.
include CMakeFiles/vulcan.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/vulcan.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vulcan.dir/flags.make

CMakeFiles/vulcan.dir/main.cpp.o: CMakeFiles/vulcan.dir/flags.make
CMakeFiles/vulcan.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vulcan.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vulcan.dir/main.cpp.o -c /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/main.cpp

CMakeFiles/vulcan.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vulcan.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/main.cpp > CMakeFiles/vulcan.dir/main.cpp.i

CMakeFiles/vulcan.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vulcan.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/main.cpp -o CMakeFiles/vulcan.dir/main.cpp.s

CMakeFiles/vulcan.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/vulcan.dir/main.cpp.o.requires

CMakeFiles/vulcan.dir/main.cpp.o.provides: CMakeFiles/vulcan.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/vulcan.dir/build.make CMakeFiles/vulcan.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/vulcan.dir/main.cpp.o.provides

CMakeFiles/vulcan.dir/main.cpp.o.provides.build: CMakeFiles/vulcan.dir/main.cpp.o


CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o: CMakeFiles/vulcan.dir/flags.make
CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o: ../src/continuum_mechanics.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o -c /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/continuum_mechanics.cpp

CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/continuum_mechanics.cpp > CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.i

CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/continuum_mechanics.cpp -o CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.s

CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.requires:

.PHONY : CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.requires

CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.provides: CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.requires
	$(MAKE) -f CMakeFiles/vulcan.dir/build.make CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.provides.build
.PHONY : CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.provides

CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.provides.build: CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o


CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o: CMakeFiles/vulcan.dir/flags.make
CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o: ../src/read_gmsh_bc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o -c /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/read_gmsh_bc.cpp

CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/read_gmsh_bc.cpp > CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.i

CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/read_gmsh_bc.cpp -o CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.s

CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.requires:

.PHONY : CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.requires

CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.provides: CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.requires
	$(MAKE) -f CMakeFiles/vulcan.dir/build.make CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.provides.build
.PHONY : CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.provides

CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.provides.build: CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o


CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o: CMakeFiles/vulcan.dir/flags.make
CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o: ../src/read_gmsh_element_property.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o -c /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/read_gmsh_element_property.cpp

CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/read_gmsh_element_property.cpp > CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.i

CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/read_gmsh_element_property.cpp -o CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.s

CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.requires:

.PHONY : CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.requires

CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.provides: CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.requires
	$(MAKE) -f CMakeFiles/vulcan.dir/build.make CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.provides.build
.PHONY : CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.provides

CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.provides.build: CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o


CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o: CMakeFiles/vulcan.dir/flags.make
CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o: ../src/vulcan_define_system.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o -c /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/vulcan_define_system.cpp

CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/vulcan_define_system.cpp > CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.i

CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/src/vulcan_define_system.cpp -o CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.s

CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.requires:

.PHONY : CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.requires

CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.provides: CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.requires
	$(MAKE) -f CMakeFiles/vulcan.dir/build.make CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.provides.build
.PHONY : CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.provides

CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.provides.build: CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o


# Object files for target vulcan
vulcan_OBJECTS = \
"CMakeFiles/vulcan.dir/main.cpp.o" \
"CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o" \
"CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o" \
"CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o" \
"CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o"

# External object files for target vulcan
vulcan_EXTERNAL_OBJECTS =

vulcan: CMakeFiles/vulcan.dir/main.cpp.o
vulcan: CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o
vulcan: CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o
vulcan: CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o
vulcan: CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o
vulcan: CMakeFiles/vulcan.dir/build.make
vulcan: /Users/gobalk/Libraries/libmesh-install/lib/libmesh_opt.dylib
vulcan: /Users/gobalk/homebrew/lib/libmpi.dylib
vulcan: /Users/gobalk/homebrew/lib/libmpi_cxx.dylib
vulcan: CMakeFiles/vulcan.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable vulcan"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vulcan.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vulcan.dir/build: vulcan

.PHONY : CMakeFiles/vulcan.dir/build

CMakeFiles/vulcan.dir/requires: CMakeFiles/vulcan.dir/main.cpp.o.requires
CMakeFiles/vulcan.dir/requires: CMakeFiles/vulcan.dir/src/continuum_mechanics.cpp.o.requires
CMakeFiles/vulcan.dir/requires: CMakeFiles/vulcan.dir/src/read_gmsh_bc.cpp.o.requires
CMakeFiles/vulcan.dir/requires: CMakeFiles/vulcan.dir/src/read_gmsh_element_property.cpp.o.requires
CMakeFiles/vulcan.dir/requires: CMakeFiles/vulcan.dir/src/vulcan_define_system.cpp.o.requires

.PHONY : CMakeFiles/vulcan.dir/requires

CMakeFiles/vulcan.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vulcan.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vulcan.dir/clean

CMakeFiles/vulcan.dir/depend:
	cd /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build /Users/gobalk/Work/vulcan/working_versions/linear_elasticity/truss/build/CMakeFiles/vulcan.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vulcan.dir/depend

