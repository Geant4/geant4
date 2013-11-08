#
# the makefile that will be used by the 'make' utility
#

#a variable which holds the name of the executable file to be created
name := radioprotection

#a variable that puts this into a variable that is used in the binmake.gmk file
#(see below)
G4TARGET := $(name)

#variable, what does it do?
G4EXLIB := true

#define the target 'all' to be 'phony'; ie, isn't a file
.PHONY: all

#target is all, dependencies are 'lib' and 'bin', ie, will check to see if anything was changed in those sub-directories
all: lib bin

 CPPFLAGS += -I${ROOTSYS}/include
 EXTRALIBS = $(shell root-config --glibs)


#now, include the "geant4" makefile, builds the source code in a standard way
include $(G4INSTALL)/config/binmake.gmk
