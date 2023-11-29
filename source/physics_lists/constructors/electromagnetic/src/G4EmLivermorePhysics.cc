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

#include "G4EmLivermorePhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models
// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPPolarizedComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreGammaConversion5DModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh"

#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"

// e+-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4Generator2BS.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4BetheHeitler5DModel.hh"

// e+
#include "G4eplusAnnihilation.hh"

// hadrons
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4ePairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
#include "G4LindhardSorensenIonModel.hh"

// msc models
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

// interfaces
#include "G4EmParameters.hh"
#include "G4EmBuilder.hh"
#include "G4EmStandUtil.hh"

// particles

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4GenericIon.hh"

//
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmLivermorePhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLivermorePhysics::G4EmLivermorePhysics(G4int ver, const G4String& pname)
  : G4VPhysicsConstructor(pname)
{
  SetVerboseLevel(ver);
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(ver);
  param->SetMinEnergy(100*CLHEP::eV);
  param->SetLowestElectronEnergy(100*CLHEP::eV);
  param->SetNumberOfBinsPerDecade(20);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetStepFunction(0.2, 10*CLHEP::um);
  param->SetStepFunctionMuHad(0.1, 50*CLHEP::um);
  param->SetStepFunctionLightIons(0.1, 20*CLHEP::um);
  param->SetStepFunctionIons(0.1, 1*CLHEP::um);
  param->SetUseMottCorrection(true);
  param->SetMscStepLimitType(fUseSafetyPlus);
  param->SetMscSkin(3);
  param->SetMscRangeFactor(0.08);
  param->SetMuHadLateralDisplacement(true);
  param->SetFluo(true);
  param->SetUseICRU90Data(true);
  param->SetFluctuationType(fUrbanFluctuation);
  param->SetMaxNIELEnergy(1*CLHEP::MeV);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLivermorePhysics::~G4EmLivermorePhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLivermorePhysics::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLivermorePhysics::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmBuilder::PrepareEMPhysics();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();
  
  // processes used by several particles
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

  // high energy limit for e+- scattering models
  G4double highEnergyLimit= param->MscEnergyLimit();
  G4double livEnergyLimit = 1*CLHEP::GeV;

  // nuclear stopping
  G4double nielEnergyLimit = param->MaxNIELEnergy();
  G4NuclearStopping* pnuc = nullptr;
  if(nielEnergyLimit > 0.0) {
    pnuc = new G4NuclearStopping();
    pnuc->SetMaxKinEnergy(nielEnergyLimit);
  }

  // Add Livermore EM Processes

  // Add gamma EM Processes
  G4ParticleDefinition* particle = G4Gamma::Gamma();
  G4bool polar = param->EnablePolarisation();

  // photoelectric effect - Livermore model only
  G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
  G4VEmModel* peModel = new G4LivermorePhotoElectricModel();
  pe->SetEmModel(peModel);
  if(polar) {
    peModel->SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized());
  }

  // Compton scattering - Livermore model only
  G4ComptonScattering* cs = new G4ComptonScattering;
  cs->SetEmModel(new G4KleinNishinaModel());
  G4VEmModel* cModel = nullptr;
  if(polar) { 
    cModel = new G4LivermorePolarizedComptonModel();
    cModel->SetHighEnergyLimit(20*CLHEP::MeV);
  } else {
    cModel = new G4LivermoreComptonModel();
    cModel->SetHighEnergyLimit(livEnergyLimit);
  }
  cs->AddEmModel(0, cModel);

  // gamma conversion 
  G4GammaConversion* gc = new G4GammaConversion();
  G4VEmModel* convLiv = new G4LivermoreGammaConversion5DModel();
  gc->SetEmModel(convLiv);

  // Rayleigh scattering
  G4RayleighScattering* rl = new G4RayleighScattering();
  if(polar) {
    rl->SetEmModel(new G4LivermorePolarizedRayleighModel());
  }

  ph->RegisterProcess(pe, particle);
  ph->RegisterProcess(cs, particle);
  ph->RegisterProcess(gc, particle);
  ph->RegisterProcess(rl, particle);

  // e-
  particle = G4Electron::Electron();

  // multiple and single scattering
  G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
  G4WentzelVIModel* msc2 = new G4WentzelVIModel();
  msc1->SetHighEnergyLimit(highEnergyLimit);
  msc2->SetLowEnergyLimit(highEnergyLimit);
  G4EmBuilder::ConstructElectronMscProcess(msc1, msc2, particle);

  G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
  G4CoulombScattering* ss = new G4CoulombScattering();
  ss->SetEmModel(ssm); 
  ss->SetMinKinEnergy(highEnergyLimit);
  ssm->SetLowEnergyLimit(highEnergyLimit);
  ssm->SetActivationLowEnergyLimit(highEnergyLimit);
      
  // Ionisation - Livermore should be used only for low energies
  G4eIonisation* eioni = new G4eIonisation();
  eioni->SetFluctModel(G4EmStandUtil::ModelOfFluctuations());
  G4VEmModel* theIoniLiv = new G4LivermoreIonisationModel();
  theIoniLiv->SetHighEnergyLimit(0.1*CLHEP::MeV); 
  eioni->AddEmModel(0, theIoniLiv);
      
  // Bremsstrahlung
  G4eBremsstrahlung* brem = new G4eBremsstrahlung();
  G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
  G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br1->SetHighEnergyLimit(GeV);

  G4ePairProduction* ee = new G4ePairProduction();
     
  // register processes
  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(ss, particle);

  // e+
  particle = G4Positron::Positron();

  // multiple and single scattering
  msc1 = new G4GoudsmitSaundersonMscModel();
  msc2 = new G4WentzelVIModel();
  msc1->SetHighEnergyLimit(highEnergyLimit);
  msc2->SetLowEnergyLimit(highEnergyLimit);
  G4EmBuilder::ConstructElectronMscProcess(msc1, msc2, particle);

  ssm = new G4eCoulombScatteringModel(); 
  ss = new G4CoulombScattering();
  ss->SetEmModel(ssm); 
  ss->SetMinKinEnergy(highEnergyLimit);
  ssm->SetLowEnergyLimit(highEnergyLimit);
  ssm->SetActivationLowEnergyLimit(highEnergyLimit);

  // ionisation
  eioni = new G4eIonisation();

  // Bremsstrahlung from standard
  brem = new G4eBremsstrahlung();
  br1 = new G4SeltzerBergerModel();
  br2 = new G4eBremsstrahlungRelModel();
  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br1->SetHighEnergyLimit(GeV);

  // register processes
  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(new G4eplusAnnihilation(), particle);
  ph->RegisterProcess(ss, particle);

  // generic ion
  particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  ionIoni->SetEmModel(new G4LindhardSorensenIonModel());
  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);
  if(nullptr != pnuc) { ph->RegisterProcess(pnuc, particle); }

  // muons, hadrons, ions
  G4EmBuilder::ConstructCharged(hmsc, pnuc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
