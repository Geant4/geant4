//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysicsWVI
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
// 13.11.2006 V.Ivanchenko use G4hMultipleScattering
// 23.11.2006 V.Ivanchenko remove mscStepLimit option and improve cout
// 13.02.2007 V.Ivanchenko use G4hMultipleScattering for muons
// 13.02.2007 V.Ivanchenko set skin=0.0
// 21.04.2008 V.Ivanchenko add long-lived D and B mesons
//
//----------------------------------------------------------------------------
//

#include "G4EmStandardPhysicsWVI.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmParameters.hh"
#include "G4EmBuilder.hh"
#include "G4LossTableManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

#include "G4KleinNishinaModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4WentzelVIModel.hh"
#include "G4WentzelVIRelModel.hh"
#include "G4hCoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4AtimaEnergyLossModel.hh"
#include "G4AtimaFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4BraggIonModel.hh"
#include "G4NuclearStopping.hh"
#include "G4eplusTo2GammaOKVIModel.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysicsWVI);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsWVI::G4EmStandardPhysicsWVI(G4int ver)
  : G4VPhysicsConstructor("G4EmStandardWVI"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  param->SetMinEnergy(10*CLHEP::eV);
  param->SetLowestElectronEnergy(10*CLHEP::eV);
  param->SetNumberOfBinsPerDecade(20);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetStepFunction(0.2, 100*CLHEP::um);
  param->SetStepFunctionMuHad(0.2, 50*CLHEP::um);
  param->SetStepFunctionLightIons(0.1, 20*CLHEP::um);
  param->SetStepFunctionIons(0.1, 1*CLHEP::um);
  param->SetUseMottCorrection(true);
  param->SetMuHadLateralDisplacement(true);
  param->SetMscThetaLimit(0.15);
  param->SetFluo(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsWVI::~G4EmStandardPhysicsWVI()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsWVI::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsWVI::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmBuilder::PrepareEMPhysics();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // common processes
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");
  G4NuclearStopping* pnuc(nullptr);

  // Add gamma EM processes
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  G4PhotoElectricEffect* pee = new G4PhotoElectricEffect();
  pee->SetEmModel(new G4LivermorePhotoElectricModel());

  G4ComptonScattering* cs = new G4ComptonScattering;
  cs->SetEmModel(new G4KleinNishinaModel());

  ph->RegisterProcess(pee, particle);
  ph->RegisterProcess(cs, particle);
  ph->RegisterProcess(new G4GammaConversion(), particle);
  ph->RegisterProcess(new G4RayleighScattering(), particle);

  // e-
  particle = G4Electron::Electron();

  G4eMultipleScattering* msc = new G4eMultipleScattering;
  msc->SetEmModel(new G4WentzelVIModel());
  G4CoulombScattering* ss = new G4CoulombScattering();

  ph->RegisterProcess(msc, particle);
  ph->RegisterProcess(new G4eIonisation(), particle);
  ph->RegisterProcess(new G4eBremsstrahlung(), particle);
  ph->RegisterProcess(ss, particle);

  // e+
  particle = G4Positron::Positron();

  msc = new G4eMultipleScattering;
  msc->SetEmModel(new G4WentzelVIModel());
  ss = new G4CoulombScattering();

  G4eplusAnnihilation* ann = new G4eplusAnnihilation();
  ann->SetEmModel(new G4eplusTo2GammaOKVIModel());

  ph->RegisterProcess(msc, particle);
  ph->RegisterProcess(new G4eIonisation(), particle);
  ph->RegisterProcess(new G4eBremsstrahlung(), particle);
  ph->RegisterProcess(ann, particle);
  ph->RegisterProcess(ss, particle);

  // generic ion
  particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  ionIoni->SetEmModel(new G4BraggIonModel(),0);
  ionIoni->SetEmModel(new G4AtimaEnergyLossModel(),1);
  ionIoni->SetFluctModel(new G4AtimaFluctuations());
  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);

  // muons, hadrons, ions
  G4EmBuilder::ConstructCharged(hmsc, pnuc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
