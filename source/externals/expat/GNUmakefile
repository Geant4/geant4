# $Id: GNUmakefile 54915 2011-11-14 17:10:43Z gcosmo $
# -------------------------------------------------------------
# Geant4 GNUmakefile for expat.  Gunter Folger 12-April-2011

name := G4expat

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifeq ($(G4SYSTEM),WIN32-VC)  # Windows
  #orig    /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_LIB" /D "COMPILED_FROM_DSP" /D "XML_UNICODE_WCHAR_T" /D "_VC80_UPGRADE=0x0600" /D "_MBCS" /FD /EHsc /RTC1 /MTd /Fp".\..\win32\tmp\Debug-w_static/expatw_static.pch"

  CPPFLAGS += -DWIN32 -D_WINDOWS -D_LIB -DCOMPILED_FROM_DSP -DXML_UNICODE_WCHAR_T -D_VC80_UPGRADE=0x0600 -D_MBCS
else                         # Linux, Mac
  CPPFLAGS += -DHAVE_EXPAT_CONFIG_H    -Isrc
endif

include $(G4INSTALL)/config/common.gmk
