# $Id: GNUmakefile,v 1.15 2003-11-10 17:50:50 trindade Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------


ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: bin

#G4EXEC_BUILD = true

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif


include $(G4INSTALL)/config/architecture.gmk

#ifdef G4ANALYSIS_USE
#  CPPFLAGS += `aida-config --include`
#  LDFLAGS  += `aida-config --lib`
#endif

#ifdef G4ANALYSIS_USE
AIDAINC = -I/afs/cern.ch/sw/contrib/AIDA/3.0/src/cpp
CPPFLAGS += $(AIDAINC)
LDFLAGS  += -L/afs/cern.ch/sw/lhcxx/specific/redhat73/gcc-2.95.2/5.0.5/lib \
-lAnaphe_AIDA_AnalysisFactory_native -lAnaphe_AIDA_Annotation_native \
-lAnaphe_AIDA_Histogram_native -lAnaphe_AIDA_Tree_native -lAnaphe_AIDA_HBookStore -lg2c-forMinuit \
-lHepUtilities -lCLHEP -lAnaphe_AIDA_Tuple_native
#endif


include $(G4INSTALL)/config/binmake.gmk




