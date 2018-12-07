# --------------------------------------------------------------
# GNUmakefile for common repository examples module. 
# --------------------------------------------------------------

name := dicomReader

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/pii/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/run/include
CPPFLAGS += -I$(G4INCLUDE)

# DCMTK support
#

#EXTRALIBS = -L$(DCMTK_BASE_DIR)/lib  -ldcmpstat -ldcmwlm -lijg8 -ldcmdata -ldcmjpeg -ldcmqrdb -li2d -loflog -ldcmdsig -ldcmjpls -ldcmsr -lijg12 -lofstd -ldcmimage -ldcmnet -ldcmtls -lijg16 -ldcmjpeg -ldcmrt -lcharls -ldcmimgle -lpthread -lpng

#CPPFLAGS += -DG4_DCMTK  -DHAVE_CONFIG_H -DUSE_NULL_SAFE_OFSTRING -DWITH_ARITHMETIC_PATCH -D_REENTRANT -D_XOPEN_SOURCE_EXTENDED -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_POSIX_C_SOURCE=199506L -fPIC
CPPFLAGS += -DG4_DCMTK

CPPFLAGS += -I$(DCMTK_BASE_DIR)/include

include $(G4INSTALL)/config/common.gmk
