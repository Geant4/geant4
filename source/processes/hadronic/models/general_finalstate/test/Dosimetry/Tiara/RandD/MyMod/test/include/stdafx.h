#ifdef G4DEBUG
#define __DEBUG__
#endif

#ifndef __NAME__
#define __MMAN_EXISTS__
#endif

#include "strings.h"
#include "ctype.h"
#include "globals.hh"
#include "signal.h"
#include "time.h"

#include "G4UImessenger.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4UserSteppingAction.hh"
#include "G4VisManager.hh"
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IVector.h"
#include "Interfaces/IHistoFactory.h"
#include "Interfaces/IVectorFactory.h"
#include "Interfaces/IPlotter.h"
#include "Interfaces/IAxis.h"

#include "cfortran.h"
#include "hplot.h"
#include "higz.h"
#include "kuip.h"

#include "AnalyzerLin.hh"
#include "NeutAnalyzer.hh"
#include "AnalyzerMessenger.hh"
#include "G4DosimPhysics.hh"
#include "GunMessenger.hh"
#include "Hall.hh"
#include "HallMessenger.hh"
#include "ParticleGun.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "VisManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ProcessManager.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "HallMessenger.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "stdlib.h"

#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4VIS_USE
#ifdef G4VIS_USE_DAWN
#include "G4FukuiRenderer.hh"
#endif
#ifdef G4VIS_USE_DAWNFILE
#include "G4DAWNFILE.hh"
#endif
#ifdef G4VIS_USE_OPACS
#include "G4Wo.hh"
#include "G4Xo.hh"
#endif
#ifdef G4VIS_USE_OPENGLX
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLStoredX.hh"
#endif
#ifdef G4VIS_USE_OPENGLWIN32
#include "G4OpenGLImmediateWin32.hh"
#include "G4OpenGLStoredWin32.hh"
#endif
#ifdef G4VIS_USE_OPENGLXM
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLStoredXm.hh"
#endif
#ifdef G4VIS_USE_OIX
#include "G4OpenInventorX.hh"
#endif
#ifdef G4VIS_USE_OIWIN32
#include "G4OpenInventorWin32.hh"
#endif
#ifdef G4VIS_USE_RAYX
#include "G4RayX.hh"
#endif
#ifdef G4VIS_USE_VRML
#include "G4VRML1.hh"
#endif
#ifdef G4VIS_USE_VRMLFILE
#include "G4VRML1File.hh"
#endif
#endif

