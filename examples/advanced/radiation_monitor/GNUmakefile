# $Id: GNUmakefile,v 1.7 2006-03-08 17:55:37 sfonseca Exp $
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

ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 CPPFLAGS += -DG4ANALYSIS_USE
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib`
endif

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
