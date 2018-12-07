# --------------------------------------------------------------
# GNUmakefile for examples module
# with use of classes from common repository.
# --------------------------------------------------------------

name := clGeometry

G4TARGET := $(name)
G4EXLIB := true
SUBDIRS := common

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all makesub clean cleansub
all: makesub lib bin

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += -I./common/include

makesub:
	@for dir in $(SUBDIRS); do ( \
	echo Entering $$dir ... ; \
	cd $$dir; \
	$(MAKE) obj name=clGeometry );\
	done

clean:: cleansub

cleansub:
	@for dir in $(SUBDIRS); do ( \
	echo Entering $$dir ...; \
	cd $$dir; \
	$(MAKE) clean );\
	done
