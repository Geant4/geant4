
name := G4RayTracer

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = libG4vis_management.lib libG4modeling.lib
GLOBLIBS += libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

# NO need OGL, then reset OGLGLAGS
OGLFLAGS =
OGLLIBS  =
# NO need QT, then reset QTGLAGS
QTFLAGS =
QTLIBS  =
GLQTLIBS =

CPPFLAGS += -DG4RAYTRACER_ALLOC_EXPORT

CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/navigation/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/geometry/magneticfield/include
CPPFLAGS += -I$(G4BASE)/digits_hits/detector/include
CPPFLAGS += -I$(G4BASE)/digits_hits/digits/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/run/include
CPPFLAGS += -I$(G4BASE)/tasking/include
CPPFLAGS += -I$(G4BASE)/track/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include
CPPFLAGS += -I$(G4BASE)/particles/bosons/include
CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/processes/management/include
CPPFLAGS += -I$(G4BASE)/processes/cuts/include
CPPFLAGS += -I$(G4BASE)/processes/scoring/include
CPPFLAGS += -I$(G4BASE)/materials/include

# Locally adjust source list before including common.gmk
sources := $(addprefix src/,G4RTJpegCoder.cc \
  G4RTJpegMaker.cc \
  G4RTMessenger.cc \
  G4RTOutBitStream.cc \
  G4RTPrimaryGeneratorAction.cc \
  G4RTRun.cc \
  G4RTRunAction.cc \
  G4RTSimpleScanner.cc \
  G4RTSteppingAction.cc \
  G4RTTrackingAction.cc \
  G4RTWorkerInitialization.cc \
  G4RayTracer.cc \
  G4RayTracerSceneHandler.cc \
  G4RayTracerViewer.cc \
  G4RayTrajectory.cc \
  G4RayTrajectoryPoint.cc \
  G4TheMTRayTracer.cc \
  G4TheRayTracer.cc \
  G4VRTScanner.cc)

ifdef G4VIS_BUILD_RAYTRACERX_DRIVER
  sources += src/G4RayTracerX.cc src/G4RayTracerXViewer.cc src/G4RTXScanner.cc
endif

include $(G4INSTALL)/config/common.gmk
