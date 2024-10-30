# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lasse/Desktop/Projects/FB_martin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lasse/Desktop/Projects/FB_martin

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/lasse/Desktop/Projects/FB_martin/CMakeFiles /home/lasse/Desktop/Projects/FB_martin//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/lasse/Desktop/Projects/FB_martin/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named MyLib

# Build rule for target.
MyLib: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 MyLib
.PHONY : MyLib

# fast build rule for target.
MyLib/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/build
.PHONY : MyLib/fast

#=============================================================================
# Target rules for targets named fwbw

# Build rule for target.
fwbw: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 fwbw
.PHONY : fwbw

# fast build rule for target.
fwbw/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fwbw.dir/build.make CMakeFiles/fwbw.dir/build
.PHONY : fwbw/fast

#=============================================================================
# Target rules for targets named cacode

# Build rule for target.
cacode: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 cacode
.PHONY : cacode

# fast build rule for target.
cacode/fast:
	$(MAKE) $(MAKESILENT) -f cacode/CMakeFiles/cacode.dir/build.make cacode/CMakeFiles/cacode.dir/build
.PHONY : cacode/fast

BaseMatrix.o: BaseMatrix.cpp.o
.PHONY : BaseMatrix.o

# target to build an object file
BaseMatrix.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/BaseMatrix.cpp.o
.PHONY : BaseMatrix.cpp.o

BaseMatrix.i: BaseMatrix.cpp.i
.PHONY : BaseMatrix.i

# target to preprocess a source file
BaseMatrix.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/BaseMatrix.cpp.i
.PHONY : BaseMatrix.cpp.i

BaseMatrix.s: BaseMatrix.cpp.s
.PHONY : BaseMatrix.s

# target to generate assembly for a file
BaseMatrix.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/BaseMatrix.cpp.s
.PHONY : BaseMatrix.cpp.s

FwBwAligner.o: FwBwAligner.cpp.o
.PHONY : FwBwAligner.o

# target to build an object file
FwBwAligner.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/FwBwAligner.cpp.o
.PHONY : FwBwAligner.cpp.o

FwBwAligner.i: FwBwAligner.cpp.i
.PHONY : FwBwAligner.i

# target to preprocess a source file
FwBwAligner.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/FwBwAligner.cpp.i
.PHONY : FwBwAligner.cpp.i

FwBwAligner.s: FwBwAligner.cpp.s
.PHONY : FwBwAligner.s

# target to generate assembly for a file
FwBwAligner.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/FwBwAligner.cpp.s
.PHONY : FwBwAligner.cpp.s

SubstitutionMatrix.o: SubstitutionMatrix.cpp.o
.PHONY : SubstitutionMatrix.o

# target to build an object file
SubstitutionMatrix.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/SubstitutionMatrix.cpp.o
.PHONY : SubstitutionMatrix.cpp.o

SubstitutionMatrix.i: SubstitutionMatrix.cpp.i
.PHONY : SubstitutionMatrix.i

# target to preprocess a source file
SubstitutionMatrix.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/SubstitutionMatrix.cpp.i
.PHONY : SubstitutionMatrix.cpp.i

SubstitutionMatrix.s: SubstitutionMatrix.cpp.s
.PHONY : SubstitutionMatrix.s

# target to generate assembly for a file
SubstitutionMatrix.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/MyLib.dir/build.make CMakeFiles/MyLib.dir/SubstitutionMatrix.cpp.s
.PHONY : SubstitutionMatrix.cpp.s

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fwbw.dir/build.make CMakeFiles/fwbw.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fwbw.dir/build.make CMakeFiles/fwbw.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/fwbw.dir/build.make CMakeFiles/fwbw.dir/main.cpp.s
.PHONY : main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... MyLib"
	@echo "... cacode"
	@echo "... fwbw"
	@echo "... BaseMatrix.o"
	@echo "... BaseMatrix.i"
	@echo "... BaseMatrix.s"
	@echo "... FwBwAligner.o"
	@echo "... FwBwAligner.i"
	@echo "... FwBwAligner.s"
	@echo "... SubstitutionMatrix.o"
	@echo "... SubstitutionMatrix.i"
	@echo "... SubstitutionMatrix.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

