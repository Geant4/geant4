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
// ClassName:   G4EmStandardPhysicsSS
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

#include "G4EmStandardPhysicsSS.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmParameters.hh"
#include "G4EmBuilder.hh"
#include "G4LossTableManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4BetheHeitler5DModel.hh"
#include "G4hMultipleScattering.hh"

#include "G4KleinNishinaModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4hCoulombScatteringModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4LindhardSorensenIonModel.hh"
#include "G4IonFluctuations.hh"

//#include "G4eSingleCoulombScatteringModel.hh"
#include "G4eDPWACoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4ePairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"
#include "G4GammaGeneralProcess.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysicsSS);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsSS::G4EmStandardPhysicsSS(G4int ver)
  : G4VPhysicsConstructor("G4EmStandardSS")
{
  SetVerboseLevel(ver);
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(ver);
  param->SetLowestElectronEnergy(10*CLHEP::eV);
  param->SetMscThetaLimit(0.0);
  param->SetUseMottCorrection(true); // use Mott-correction for e-/e+ msc gs
  param->SetAuger(true);
  param->SetPixe(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysicsSS::~G4EmStandardPhysicsSS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsSS::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysicsSS::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmBuilder::PrepareEMPhysics();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();

  // processes used by several particles
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

  // Add gamma EM Processes
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  // Photoelectric
  G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
  G4VEmModel* peModel = new G4LivermorePhotoElectricModel();
  pe->SetEmModel(peModel);
  if(param->EnablePolarisation()) {
    peModel->SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized());
  }

  // Compton scattering
  G4ComptonScattering* cs = new G4ComptonScattering;
  cs->SetEmModel(new G4KleinNishinaModel());

  // Gamma conversion
  G4GammaConversion* gc = new G4GammaConversion();
  G4VEmModel* conv = new G4BetheHeitler5DModel();
  gc->SetEmModel(conv);

  // default Rayleigh scattering is Livermore
  G4RayleighScattering* rl = new G4RayleighScattering();
  if(param->EnablePolarisation()) {
    rl->SetEmModel(new G4LivermorePolarizedRayleighModel());
  }

  if(param->GeneralProcessActive()) {
    G4GammaGeneralProcess* sp = new G4GammaGeneralProcess();
    sp->AddEmProcess(pe);
    sp->AddEmProcess(cs);
    sp->AddEmProcess(gc);
    sp->AddEmProcess(rl);
    G4LossTableManager::Instance()->SetGammaGeneralProcess(sp);
    ph->RegisterProcess(sp, particle);
  } else {
    ph->RegisterProcess(pe, particle);
    ph->RegisterProcess(cs, particle);
    ph->RegisterProcess(gc, particle);
    ph->RegisterProcess(rl, particle);
  }
  // e-
  particle = G4Electron::Electron();

  G4CoulombScattering* ss = new G4CoulombScattering();
  if(param->UseMottCorrection()) {
    ss->SetEmModel(new G4eDPWACoulombScatteringModel());
  } else {
    ss->SetEmModel(new G4eCoulombScatteringModel()); 
  }
  ph->RegisterProcess(new G4eIonisation(), particle);
  ph->RegisterProcess(new G4eBremsstrahlung(), particle);

  G4ePairProduction* ee = new G4ePairProduction();
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(ss, particle);

  // e+
  particle = G4Positron::Positron();

  ss = new G4CoulombScattering();
  if(param->UseMottCorrection()) {
    ss->SetEmModel(new G4eDPWACoulombScatteringModel());
  } else {
    ss->SetEmModel(new G4eCoulombScatteringModel()); 
  }
  ph->RegisterProcess(new G4eIonisation(), particle);
  ph->RegisterProcess(new G4eBremsstrahlung(), particle);
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(ss, particle);
  ph->RegisterProcess(new G4eplusAnnihilation(), particle);

  // generic ion
  particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  auto fluc = new G4IonFluctuations();
  ionIoni->SetFluctModel(fluc);
  ionIoni->SetEmModel(new G4LindhardSorensenIonModel());
  ph->RegisterProcess(ionIoni, particle);
  ph->RegisterProcess(new G4CoulombScattering(), particle);

  // muons, hadrons, ions
  G4EmBuilder::ConstructChargedSS(hmsc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
