# $Id: GNUmakefile,v 1.8 2009-03-06 09:32:29 gbarrand Exp $
# --------------------------------------------------------------
# GNUmakefile for examples with internal subpackages
# R. Capra, 09 Sep 2005
# --------------------------------------------------------------

name := Radmon
G4TARGET := $(name)
G4EXLIB := true
SUBDIRS := userActions generalPurpose materials sensitivity detector physics generator 

ifdef G4ANALYSIS_USE
 SUBDIRS += analysis
endif

SUBDIRS += application

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all bin makesub clean cleansub
all: makesub lib bin

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += -I./detector/include \
            -I./generalPurpose/include \
            -I./application/include

makesub:
	@for dir in $(SUBDIRS); do ( \
		echo Entering $$dir ... ; \
		cd $$dir; \
		$(MAKE) obj );\
	done

clean:: cleansub

cleansub:
	@for dir in $(SUBDIRS); do ( \
		echo Entering $$dir ...; \
		cd $$dir; \
		$(MAKE) clean );\
	done
