# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := DICOM
G4TARGET := $(name)
G4EXLIB := true
SUBDIRS := dicomReader

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all  makesub clean cleansub
all:  makesub lib bin

# DCMTK support
#
#DICOM_USE_DCMTK := true

ifdef DICOM_USE_DCMTK

EXTRALIBS = -L$(DCMTK_BASE_DIR)/lib  -ldcmpstat -ldcmwlm -lijg8 -ldcmdata -ldcmjpeg -ldcmqrdb -li2d -loflog -ldcmdsig -ldcmjpls -ldcmsr -lijg12 -lofstd -ldcmimage -ldcmnet -ldcmtls -lijg16 -ldcmjpeg -ldcmrt -lcharls -ldcmimgle -lpthread -lpng

#CPPFLAGS += -DG4_DCMTK  -DHAVE_CONFIG_H -DUSE_NULL_SAFE_OFSTRING -DWITH_ARITHMETIC_PATCH -D_REENTRANT -D_XOPEN_SOURCE_EXTENDED -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_POSIX_C_SOURCE=199506L -fPIC
CPPFLAGS += -DG4_DCMTK 

endif

include $(G4INSTALL)/config/binmake.gmk

ifdef DICOM_USE_DCMTK
INCFLAGS += -I$(DCMTK_BASE_DIR)/include
CPPFLAGS += -I./dicomReader/include
endif

makesub:
ifdef DICOM_USE_DCMTK
	@for dir in $(SUBDIRS); do ( \
	echo Entering $$dir ... ; \
	cd $$dir; \
	$(MAKE) obj name=dicomReader );\
	done
endif

clean:: cleansub

cleansub:
ifdef DICOM_USE_DCMTK
	@for dir in $(SUBDIRS); do ( \
	echo Entering $$dir ...; \
	cd $$dir; \
	$(MAKE) clean );\
	done
endif

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
