# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sakura/Documents/CS174/graphicsMatte

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sakura/Documents/CS174/graphicsMatte

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/sakura/Documents/CS174/graphicsMatte/CMakeFiles /home/sakura/Documents/CS174/graphicsMatte/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/sakura/Documents/CS174/graphicsMatte/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named figtree

# Build rule for target.
figtree: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 figtree
.PHONY : figtree

# fast build rule for target.
figtree/fast:
	$(MAKE) -f CMakeFiles/figtree.dir/build.make CMakeFiles/figtree.dir/build
.PHONY : figtree/fast

#=============================================================================
# Target rules for targets named geodesics

# Build rule for target.
geodesics: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 geodesics
.PHONY : geodesics

# fast build rule for target.
geodesics/fast:
	$(MAKE) -f CMakeFiles/geodesics.dir/build.make CMakeFiles/geodesics.dir/build
.PHONY : geodesics/fast

#=============================================================================
# Target rules for targets named kcc

# Build rule for target.
kcc: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 kcc
.PHONY : kcc

# fast build rule for target.
kcc/fast:
	$(MAKE) -f CMakeFiles/kcc.dir/build.make CMakeFiles/kcc.dir/build
.PHONY : kcc/fast

#=============================================================================
# Target rules for targets named kde

# Build rule for target.
kde: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 kde
.PHONY : kde

# fast build rule for target.
kde/fast:
	$(MAKE) -f CMakeFiles/kde.dir/build.make CMakeFiles/kde.dir/build
.PHONY : kde/fast

#=============================================================================
# Target rules for targets named matting

# Build rule for target.
matting: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 matting
.PHONY : matting

# fast build rule for target.
matting/fast:
	$(MAKE) -f CMakeFiles/matting.dir/build.make CMakeFiles/matting.dir/build
.PHONY : matting/fast

matting.o: matting.cpp.o
.PHONY : matting.o

# target to build an object file
matting.cpp.o:
	$(MAKE) -f CMakeFiles/matting.dir/build.make CMakeFiles/matting.dir/matting.cpp.o
.PHONY : matting.cpp.o

matting.i: matting.cpp.i
.PHONY : matting.i

# target to preprocess a source file
matting.cpp.i:
	$(MAKE) -f CMakeFiles/matting.dir/build.make CMakeFiles/matting.dir/matting.cpp.i
.PHONY : matting.cpp.i

matting.s: matting.cpp.s
.PHONY : matting.s

# target to generate assembly for a file
matting.cpp.s:
	$(MAKE) -f CMakeFiles/matting.dir/build.make CMakeFiles/matting.dir/matting.cpp.s
.PHONY : matting.cpp.s

src/KCenterClustering.o: src/KCenterClustering.cpp.o
.PHONY : src/KCenterClustering.o

# target to build an object file
src/KCenterClustering.cpp.o:
	$(MAKE) -f CMakeFiles/kcc.dir/build.make CMakeFiles/kcc.dir/src/KCenterClustering.cpp.o
.PHONY : src/KCenterClustering.cpp.o

src/KCenterClustering.i: src/KCenterClustering.cpp.i
.PHONY : src/KCenterClustering.i

# target to preprocess a source file
src/KCenterClustering.cpp.i:
	$(MAKE) -f CMakeFiles/kcc.dir/build.make CMakeFiles/kcc.dir/src/KCenterClustering.cpp.i
.PHONY : src/KCenterClustering.cpp.i

src/KCenterClustering.s: src/KCenterClustering.cpp.s
.PHONY : src/KCenterClustering.s

# target to generate assembly for a file
src/KCenterClustering.cpp.s:
	$(MAKE) -f CMakeFiles/kcc.dir/build.make CMakeFiles/kcc.dir/src/KCenterClustering.cpp.s
.PHONY : src/KCenterClustering.cpp.s

src/figtree.o: src/figtree.cpp.o
.PHONY : src/figtree.o

# target to build an object file
src/figtree.cpp.o:
	$(MAKE) -f CMakeFiles/figtree.dir/build.make CMakeFiles/figtree.dir/src/figtree.cpp.o
.PHONY : src/figtree.cpp.o

src/figtree.i: src/figtree.cpp.i
.PHONY : src/figtree.i

# target to preprocess a source file
src/figtree.cpp.i:
	$(MAKE) -f CMakeFiles/figtree.dir/build.make CMakeFiles/figtree.dir/src/figtree.cpp.i
.PHONY : src/figtree.cpp.i

src/figtree.s: src/figtree.cpp.s
.PHONY : src/figtree.s

# target to generate assembly for a file
src/figtree.cpp.s:
	$(MAKE) -f CMakeFiles/figtree.dir/build.make CMakeFiles/figtree.dir/src/figtree.cpp.s
.PHONY : src/figtree.cpp.s

src/geodesics.o: src/geodesics.cpp.o
.PHONY : src/geodesics.o

# target to build an object file
src/geodesics.cpp.o:
	$(MAKE) -f CMakeFiles/geodesics.dir/build.make CMakeFiles/geodesics.dir/src/geodesics.cpp.o
.PHONY : src/geodesics.cpp.o

src/geodesics.i: src/geodesics.cpp.i
.PHONY : src/geodesics.i

# target to preprocess a source file
src/geodesics.cpp.i:
	$(MAKE) -f CMakeFiles/geodesics.dir/build.make CMakeFiles/geodesics.dir/src/geodesics.cpp.i
.PHONY : src/geodesics.cpp.i

src/geodesics.s: src/geodesics.cpp.s
.PHONY : src/geodesics.s

# target to generate assembly for a file
src/geodesics.cpp.s:
	$(MAKE) -f CMakeFiles/geodesics.dir/build.make CMakeFiles/geodesics.dir/src/geodesics.cpp.s
.PHONY : src/geodesics.cpp.s

src/kde.o: src/kde.cpp.o
.PHONY : src/kde.o

# target to build an object file
src/kde.cpp.o:
	$(MAKE) -f CMakeFiles/kde.dir/build.make CMakeFiles/kde.dir/src/kde.cpp.o
.PHONY : src/kde.cpp.o

src/kde.i: src/kde.cpp.i
.PHONY : src/kde.i

# target to preprocess a source file
src/kde.cpp.i:
	$(MAKE) -f CMakeFiles/kde.dir/build.make CMakeFiles/kde.dir/src/kde.cpp.i
.PHONY : src/kde.cpp.i

src/kde.s: src/kde.cpp.s
.PHONY : src/kde.s

# target to generate assembly for a file
src/kde.cpp.s:
	$(MAKE) -f CMakeFiles/kde.dir/build.make CMakeFiles/kde.dir/src/kde.cpp.s
.PHONY : src/kde.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... figtree"
	@echo "... geodesics"
	@echo "... kcc"
	@echo "... kde"
	@echo "... matting"
	@echo "... rebuild_cache"
	@echo "... matting.o"
	@echo "... matting.i"
	@echo "... matting.s"
	@echo "... src/KCenterClustering.o"
	@echo "... src/KCenterClustering.i"
	@echo "... src/KCenterClustering.s"
	@echo "... src/figtree.o"
	@echo "... src/figtree.i"
	@echo "... src/figtree.s"
	@echo "... src/geodesics.o"
	@echo "... src/geodesics.i"
	@echo "... src/geodesics.s"
	@echo "... src/kde.o"
	@echo "... src/kde.i"
	@echo "... src/kde.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

