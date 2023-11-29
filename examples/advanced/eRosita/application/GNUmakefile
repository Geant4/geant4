# --------------------------------------------------------------------
# GNUmakefile for example subpackage
# --------------------------------------------------------------------

name := eRosita

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_USE.gmk
include $(G4INSTALL)/config/G4VIS_USE.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
              -I$(G4BASE)/global/HEPRandom/include \
              -I$(G4BASE)/global/HEPGeometry/include \
              -I$(G4BASE)/global/HEPNumerics/include \
              -I$(G4BASE)/run/include \
              -I$(G4BASE)/materials/include \
              -I$(G4BASE)/parameterisations/trd_clusters/include \
              -I$(G4BASE)/parameterisations/gflash/include \
              -I$(G4BASE)/particles/management/include \
              -I$(G4BASE)/particles/bosons/include \
              -I$(G4BASE)/particles/leptons/include \
              -I$(G4BASE)/particles/hadrons/barions/include \
              -I$(G4BASE)/particles/hadrons/ions/include \
              -I$(G4BASE)/particles/hadrons/mesons/include \
              -I$(G4BASE)/particles/shortlived/include \
              -I$(G4BASE)/processes/management/include \
              -I$(G4BASE)/processes/cuts/include \
              -I$(G4BASE)/processes/decay/include \
              -I$(G4BASE)/processes/optical/include \
              -I$(G4BASE)/processes/transportation/include \
              -I$(G4BASE)/processes/parameterisation/include \
              -I$(G4BASE)/processes/electromagnetic/standard/include \
              -I$(G4BASE)/processes/electromagnetic/muons/include \
              -I$(G4BASE)/processes/electromagnetic/utils/include \
              -I$(G4BASE)/processes/electromagnetic/xrays/include \
              -I$(G4BASE)/processes/electromagnetic/highenergy/include \
              -I$(G4BASE)/processes/electromagnetic/lowenergy/include \
              -I$(G4BASE)/processes/electromagnetic/pii/include \
	      -I../physics/include \
              -I$(G4BASE)/geometry/management/include \
              -I$(G4BASE)/geometry/solids/CSG/include \
              -I$(G4BASE)/geometry/solids/specific/include \
              -I$(G4BASE)/geometry/solids/BREPS/include \
              -I$(G4BASE)/geometry/solids/Boolean/include \
              -I$(G4BASE)/geometry/divisions/include \
              -I$(G4BASE)/geometry/volumes/include \
              -I$(G4BASE)/geometry/navigation/include \
              -I$(G4BASE)/geometry/magneticfield/include \
              -I$(G4BASE)/geometry/biasing/include \
              -I$(G4BASE)/track/include \
              -I$(G4BASE)/tracking/include \
              -I$(G4BASE)/digits_hits/detector/include \
              -I$(G4BASE)/digits_hits/hits/include \
              -I$(G4BASE)/digits_hits/digits/include \
              -I$(G4BASE)/digits_hits/utils/include \
              -I$(G4BASE)/persistency/include \
              -I$(G4BASE)/readout/include \
              -I$(G4BASE)/event/include \
              -I$(G4BASE)/graphics_reps/include \
              -I$(G4BASE)/intercoms/include \

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
  LOADLIBS += `aida-config --lib`
endif

include $(G4INSTALL)/config/common.gmk
