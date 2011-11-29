# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.2 2010-11-29 10:00:05 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------


name := eRosita

G4TARGET := $(name)
G4EXLIB := true
SUBDIRS := physics application

# ifdef G4ANALYSIS_USE
#  SUBDIRS += analysis
# endif

ifndef G4INSTALL
  G4INSTALL = ../..
endif

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib`
endif

CPPFLAGS += -I./physics/include \
            -I./application/include
            
.PHONY: all makesub clean cleansub
all: makesub lib bin

include $(G4INSTALL)/config/binmake.gmk

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
