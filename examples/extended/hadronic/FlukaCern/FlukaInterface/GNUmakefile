########################################
########## HOW TO USE
########################################

######### Clean build
# make clean

######### Build evth: both interface to G4 and its tests:
# make
# (or make all)
# (obviously, you can also do: make -j8)

######### Build test executables only
# make test

######### Build interface to G4 only
# make interface

######### Create environment scripts
# make env
# (creates env_FLUKA.sh and env_FLUKA_G4_interface.sh)

######### Create env_FLUKA_G4_interface.sh only
# make env_FLUKA_G4_interface
# (env_FLUKA_G4_interface.sh needs to be sourced before compiling / running any test or G4 application using the FLUKA-G4 interface).



########################################
########## COMPILER & FLAGS
########################################

######### C++ compiler and flags
CXX_COMP := g++
CXX_FLAGS := -std=c++17 -O2 -g -Wall -pedantic -malign-double


######### Fortran compiler and flags
F_COMP := gfortran
F_FLAGS := -g -Wall -pedantic


######### Linker flag
LD_FLAGS := -lgfortran


######### Check compilers collection version
IS_CXX_COMP_GT_7 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 7)
ifneq "$(IS_CXX_COMP_GT_7)" "1"
$(error "gcc suite version = $(shell gcc -dumpversion). Only versions >= 7 are supported!")
endif



########################################
########## EXTERNAL DEPENDENCIES
########################################

######### FLUKA
FLUKA_PATH := $(realpath $(shell dirname $(shell which fluka))/..)
ifeq ($(wildcard $(FLUKA_PATH)/.*),)
$(error "fluka executables not found. FLUKA_REPO=$(FLUKA_PATH). Add your fluka_repo/bin to the PATH." )
endif

FLUKA_INCLUDE_DIR := $(FLUKA_PATH)/include
ifeq ($(wildcard $(FLUKA_INCLUDE_DIR)/.*),)
$(error "FLUKA_INCLUDE_DIR=$(FLUKA_INCLUDE_DIR)/ does not contain any file. Add your fluka_repo/bin to the PATH." )
endif
FLUKA_LIB_DIR := $(FLUKA_PATH)/lib
ifeq ($(wildcard $(FLUKA_LIB_DIR)/$(FLUKA_LIB)),)
$(error "FLUKA_LIB_DIR=$(FLUKA_LIB_DIR)/ does not contain $(FLUKA_LIB). Add your fluka_repo/bin to the PATH." )
endif
FLUKA_LIB_NAME := fluka
FLUKA_LIB := lib$(FLUKA_LIB_NAME).a
FLUKA_NO_RANDOM_LIB_NAME := fluka_no_random
FLUKA_NO_RANDOM_LIB := lib$(FLUKA_NO_RANDOM_LIB_NAME).a
FLUKA_OBJ_REMOVE := flrm64.o flrndm.o flrnlp.o flrnoc.o
FLUKA_RELEASE_LIB_LINK := -L$(FLUKA_LIB_DIR) -l$(FLUKA_LIB_NAME)
FLUKA_LIB_LINK := -L$(FLUKA_LIB_DIR) -l$(FLUKA_NO_RANDOM_LIB_NAME)


######### G4
G4_PATH := $(realpath $(shell geant4-config --prefix))
G4_INCLUDE_DIR := $(G4_PATH)/include/Geant4
G4_INCLUDE := -I$(G4_INCLUDE_DIR)
G4_CXX_FLAGS := $(shell geant4-config --cflags)
G4_LIBRARIES := $(shell geant4-config --libs)



########################################
########## PROJECT DIRECTORIES
########################################

# Base dir
BASE_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

# Source code
SRC_DIR := $(BASE_DIR)/source

# Directory where to put object and dependency files
OBJ_DIR := $(BASE_DIR)/obj
INTERFACE_OBJ_DIR := $(OBJ_DIR)/interface
TEST_INTERFACE_TO_G4_OBJ_DIR := $(OBJ_DIR)/test_interface_to_G4

# Directory where to put executables
BIN_DIR := $(BASE_DIR)/bin

CLEAN_DIR := $(OBJ_DIR) $(BIN_DIR)



########################################
########## PROJECT FILES
########################################

######### FORTRAN SOURCE
F_SRC_DIR := $(SRC_DIR)/fluka4_wrapper/fortran/procedures \
	$(SRC_DIR)/cpp_fortran_bridges/fortran
	
F_SRC := $(foreach src_dir, $(F_SRC_DIR), $(wildcard $(src_dir)/*.f))
F_OBJ := $(addprefix $(OBJ_DIR)/, $(notdir $(F_SRC:.f=.o)))

######### C++ SOURCE
CXX_SRC_DIR := $(SRC_DIR)/fluka4_wrapper/cpp/procedures \
	$(SRC_DIR)/cpp_fortran_bridges/cpp \
	$(SRC_DIR)/cpp_utils
	
CXX_INCLUDE_DIR := $(CXX_SRC_DIR) $(FLUKA_INCLUDE_DIR)
CXX_INCLUDE := $(addprefix -I, $(CXX_INCLUDE_DIR))
CXX_SRC := $(foreach src_dir, $(CXX_SRC_DIR), $(wildcard $(src_dir)/*.cc))
CXX_OBJ := $(addprefix $(OBJ_DIR)/, $(notdir $(CXX_SRC:.cc=.o)))

######### INTERFACE TO G4
INTERFACE_SRC_DIR := $(SRC_DIR)/fluka5/hadronic_interactions/interface_to_G4 \
	$(SRC_DIR)/fluka_G4_bridges \
	$(SRC_DIR)/fluka_G4_physics_list
	
INTERFACE_INCLUDE_DIR := $(INTERFACE_SRC_DIR)
INTERFACE_INCLUDE := $(addprefix -I, $(INTERFACE_INCLUDE_DIR))
INTERFACE_SRC := $(foreach src_dir, $(INTERFACE_SRC_DIR), $(wildcard $(src_dir)/*.cc))
INTERFACE_OBJ := $(addprefix $(INTERFACE_OBJ_DIR)/, $(notdir $(INTERFACE_SRC:.cc=.o)))

######### TEST INTERFACE TO G4 (EXECUTABLE)
TEST_INTERFACE_TO_G4_SRC_DIR := $(SRC_DIR)/fluka5/hadronic_interactions/test/interface_to_G4
TEST_INTERFACE_TO_G4_SRC := $(foreach src_dir, $(TEST_INTERFACE_TO_G4_SRC_DIR), $(wildcard $(src_dir)/*.cc))
TEST_INTERFACE_TO_G4_INCLUDE_DIR := $(TEST_INTERFACE_TO_G4_SRC_DIR)
TEST_INTERFACE_TO_G4_INCLUDE := $(addprefix -I, $(TEST_INTERFACE_TO_G4_INCLUDE_DIR))
TEST_INTERFACE_TO_G4_OBJ := $(addprefix $(TEST_INTERFACE_TO_G4_OBJ_DIR)/, $(notdir $(TEST_INTERFACE_TO_G4_SRC:.cc=.o)))
TEST_INTERFACE_TO_G4_EXE := $(BIN_DIR)/test_interface



########################################
########## TARGETS RULES
########################################

.PHONY: all
all: interface test
	@echo
	@echo "Full build successfull."
	
.PHONY: test
test: test_interface_to_g4
	@echo
	@echo "Full tests build successfull."
	
.PHONY: interface
interface: check_G4 $(INTERFACE_OBJ) fluka_no_random_lib
	@echo
	@echo "Full FLUKA-G4 interface build successfull."
	
.PHONY: check_G4
check_G4:
	@[ -d $(G4_PATH) ] || { echo "!!! G4 repo not found. Please source geant4.sh, or G4 setup.sh (need to have G4 already installed). !!!"; exit 1; }
	@[ -d $(G4_INCLUDE_DIR) ] || { echo "!!! G4_INCLUDE_DIR=$(G4_INCLUDE_DIR)/ does not exist. !!! "; exit 1; }
	
.PHONY: env
env: env_FLUKA env_FLUKA_G4_interface
	@echo
	@echo "Full env scripts creation successfull."
	
.PHONY: env_FLUKA
env_FLUKA: env_FLUKA.sh
env_FLUKA.sh: GNUmakefile
	@echo '#! /bin/bash'                                                    >  $@
	@echo                                                                   >> $@
	@echo 'export FLUKA_PATH=$$(realpath $$(dirname $$(which fluka))/..)'   >> $@
	@echo 'if [ ! -d "$$FLUKA_PATH" ]; then'                                >> $@
	@echo '   echo "!!! fluka not found. FLUKA_REPO=$$FLUKA_PATH. Add your fluka_repo/bin to the PATH. !!!";'  >> $@
	@echo 'fi'                                                              >> $@
	@echo                                                                   >> $@
	@echo 'export FLUKADATA=$$FLUKA_PATH/data'                              >> $@
	@echo 'if [ ! -d "$$FLUKADATA" ]; then'                                 >> $@
	@echo '   echo "!!! fluka data repo not found. FLUKADATA=$$FLUKADATA. Add your fluka_repo/bin to the PATH. !!!";'  >> $@
	@echo 'fi'                                                              >> $@
	@echo                                                                   >> $@
	@echo 'export FLUKA_LIB_DIR=$$FLUKA_PATH/lib'                           >> $@
	@echo 'if [ ! -d "$$FLUKA_LIB_DIR" ]; then'                             >> $@
	@echo '   echo "!!! fluka lib repo not found. FLUKA_LIB_DIR=$$FLUKA_LIB_DIR. Add your fluka_repo/bin to the PATH. !!!";'  >> $@
	@echo 'fi'                                                              >> $@
	@echo
	@echo "Created file: env_FLUKA.sh"
	
.PHONY: env_FLUKA_G4_interface
env_FLUKA_G4_interface: env_FLUKA_G4_interface.sh
env_FLUKA_G4_interface.sh: GNUmakefile
	@echo '#! /bin/bash'                                                    >  $@
	@echo 'echo "source $(BASE_DIR)/env_FLUKA.sh"'                          >> $@
	@echo 'if [ ! -e "$(BASE_DIR)/env_FLUKA.sh" ]; then'                    >> $@
	@echo '   echo "env_FLUKA.sh does not exist. Please create it: go to $(BASE_DIR), and do make env. Then, source again this script."'   >> $@
	@echo 'else'                                                            >> $@
	@echo '   source $(BASE_DIR)/env_FLUKA.sh'                              >> $@
	@echo '   echo ""'                                                      >> $@
	@echo '   echo "Setting FLUKAInterface_INCLUDE_DIR:"'                  >> $@
	@echo '   export FLUKAInterface_INCLUDE_DIR="$(SRC_DIR)/cpp_utils $(INTERFACE_INCLUDE_DIR)"'  >> $@
	@echo '   echo "export FLUKAInterface_INCLUDE_DIR=$$FLUKAInterface_INCLUDE_DIR"'  >> $@
	@echo '   echo "FLUKAInterface_INCLUDE_DIR successfully set."'         >> $@
	@echo '   echo ""'                                                      >> $@
	@echo '   echo "Setting FLUKAInterface_LIBRARIES:"'                    >> $@
	@echo '   export FLUKAInterface_LIBRARIES="$(INTERFACE_OBJ) $(CXX_OBJ) $(F_OBJ)"'  >> $@
	@echo '   missing_obj=false;'                                           >> $@
	@echo '   for obj in $$FLUKAInterface_LIBRARIES; do'                   >> $@
	@echo '      if [ ! -e "$$obj" ]; then'                                 >> $@
	@echo '         echo "$$obj does not exist."'                           >> $@
	@echo '         missing_obj=true'                                       >> $@
	@echo '      fi'                                                        >> $@
	@echo '   done'                                                         >> $@
	@echo '   if [ "$$missing_obj" = true ]; then'                          >> $@
	@echo '      echo "!!! Object file(s) in FLUKAInterface_LIBRARIES is/are missing. Please (re)compile the G4 interface (make or make interface)!";'  >> $@
	@echo '   else'                                                         >> $@
	@echo '      export FLUKAInterface_LIBRARIES="$$FLUKAInterface_LIBRARIES $(FLUKA_LIB_DIR)/$(FLUKA_NO_RANDOM_LIB)"'  >> $@
	@echo '      export FLUKAInterface_LIBRARIES="$$FLUKAInterface_LIBRARIES $(LD_FLAGS)"'  >> $@
	@echo '      echo "export FLUKAInterface_LIBRARIES=$$FLUKAInterface_LIBRARIES"'  >> $@
	@echo '      echo "FLUKAInterface_LIBRARIES successfully set."'        >> $@
	@echo '   fi'                                                           >> $@
	@echo 'fi'                                                              >> $@
	@echo 'echo ""'                                                         >> $@
	@echo "Created file: env_FLUKA_G4_interface.sh"

.PHONY: clean
clean:
	rm -fr $(CLEAN_DIR) env_FLUKA_G4_interface.sh



########################################
########## COMPILE & LINK
########################################

# BUILD TEST INTERFACE TO G4 (EXECUTABLE)
.PHONY: test_interface_to_g4
test_interface_to_g4: fluka_no_random_lib $(TEST_INTERFACE_TO_G4_EXE)
$(TEST_INTERFACE_TO_G4_EXE): $(TEST_INTERFACE_TO_G4_OBJ) $(INTERFACE_OBJ) $(CXX_OBJ) $(F_OBJ)
	@[ -d $(@D) ] || mkdir -p $(@D)
	@echo
	@echo "Building $(notdir $@) executable to test FLUKA-G4 interface:"
	$(CXX_COMP) $^ $(FLUKA_LIB_LINK) $(LD_FLAGS) $(G4_LIBRARIES) -o $@
	@echo "Built $(notdir $@) successfully."
	
# COMPILE TEST INTERFACE TO G4 SOURCES
define TEST_INTERFACE_TO_G4_BUILD_template
$(TEST_INTERFACE_TO_G4_OBJ_DIR)/%.o: $(1)/%.cc
	@echo
	@echo "Building $$(notdir $$@) to test FLUKA-G4 interface:"
	@[ -d $$(@D) ] || mkdir -p $$(@D)
	$(CXX_COMP) -DG4_USE_FLUKA $(G4_CXX_FLAGS) $(CXX_INCLUDE) $(TEST_INTERFACE_TO_G4_INCLUDE) $(INTERFACE_INCLUDE) -c $$< -o $$@
	@echo "Built $$(notdir $$@) successfully."
endef
$(foreach src_dir, $(TEST_INTERFACE_TO_G4_SRC_DIR), $(eval $(call TEST_INTERFACE_TO_G4_BUILD_template, $(src_dir))))

# COMPILE INTERFACE SOURCES
define INTERFACE_BUILD_template
$(INTERFACE_OBJ_DIR)/%.o: $(1)/%.cc $(CXX_OBJ) $(F_OBJ)
	@echo
	@echo "Building $$(notdir $$@) for FLUKA-G4 interface:"
	@[ -d $$(@D) ] || mkdir -p $$(@D)
	$(CXX_COMP) -DG4_USE_FLUKA $(G4_CXX_FLAGS) $(CXX_INCLUDE) $(INTERFACE_INCLUDE) -c $$< -o $$@
	@echo "Built $$(notdir $$@) successfully."
endef
$(foreach src_dir, $(INTERFACE_SRC_DIR), $(eval $(call INTERFACE_BUILD_template, $(src_dir))))

# COMPILE CPP SOURCE
define CXX_BUILD_template
$(OBJ_DIR)/%.o: $(1)/%.cc
	@echo
	@echo "Building $$(notdir $$@) for core FLUKA5 C++ source:"
	@[ -d $$(@D) ] || mkdir -p $$(@D)
	$(CXX_COMP) -DG4_USE_FLUKA $(CXX_FLAGS) -c $$< -o $$@
	@echo "Built $$(notdir $$@) successfully."
endef
$(foreach src_dir, $(CXX_SRC_DIR), $(eval $(call CXX_BUILD_template, $(src_dir))))

# COMPILE FORTRAN SOURCE
define F_BUILD_template
$(OBJ_DIR)/%.o: $(1)/%.f
	@echo
	@echo "Building $$(notdir $$@) for core FLUKA5 FORTRAN source:"
	@[ -d $$(@D) ] || mkdir -p $$(@D)
	$(F_COMP) -DG4_USE_FLUKA $(F_FLAGS) -c $$< -o $$@
	@echo "Built $$(notdir $$@) successfully."
endef
$(foreach src_dir, $(F_SRC_DIR), $(eval $(call F_BUILD_template, $(src_dir))))

	

########################################
########## FLUKA4 WAS UPDATED	
########################################

.PHONY: fluka
fluka: fluka_no_random_lib

.PHONY: fluka_no_random_lib
fluka_no_random_lib: $(FLUKA_LIB_DIR)/$(FLUKA_NO_RANDOM_LIB)
$(FLUKA_LIB_DIR)/$(FLUKA_NO_RANDOM_LIB): $(FLUKA_LIB_DIR)/$(FLUKA_LIB)
	@echo
	@cp $(FLUKA_LIB_DIR)/$(FLUKA_LIB) $(FLUKA_LIB_DIR)/$(FLUKA_NO_RANDOM_LIB)
	@echo "cp $(FLUKA_LIB_DIR)/$(FLUKA_LIB) $(FLUKA_LIB_DIR)/$(FLUKA_NO_RANDOM_LIB)"
	$(foreach obj,$(FLUKA_OBJ_REMOVE),ar d $(FLUKA_LIB_DIR)/$(FLUKA_NO_RANDOM_LIB) $(obj);)
	@echo "$(FLUKA_NO_RANDOM_LIB) creation / update successfull."

