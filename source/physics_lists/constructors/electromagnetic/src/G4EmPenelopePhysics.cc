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

#include "G4EmPenelopePhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// *** Processes and models

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4PenelopePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4PenelopeComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4PenelopeGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4PenelopeRayleighModel.hh"

#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"

// e- and e+
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4PenelopeIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4PenelopeBremsstrahlungModel.hh"

// e+ only
#include "G4eplusAnnihilation.hh"
#include "G4PenelopeAnnihilationModel.hh"

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
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

// utilities
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
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmPenelopePhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmPenelopePhysics::G4EmPenelopePhysics(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmPenelope")
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
  param->SetPIXEElectronCrossSectionModel("Penelope");
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmPenelopePhysics::~G4EmPenelopePhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmPenelopePhysics::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmPenelopePhysics::ConstructProcess()
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
  G4double highEnergyLimit = param->MscEnergyLimit();

  // nuclear stopping
  G4double nielEnergyLimit = param->MaxNIELEnergy();
  G4NuclearStopping* pnuc = nullptr;
  if(nielEnergyLimit > 0.0) {
    pnuc = new G4NuclearStopping();
    pnuc->SetMaxKinEnergy(nielEnergyLimit);
  }

  //Applicability range for Penelope models
  //for higher energies, the Standard models are used   
  G4double PenelopeHighEnergyLimit = 1.0*CLHEP::GeV;

  // Add gamma EM Processes
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  //Photo-electric effect
  G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
  thePhotoElectricEffect->SetEmModel(new G4PEEffectFluoModel());
  G4PenelopePhotoElectricModel* thePEPenelopeModel = new G4PenelopePhotoElectricModel();   
  thePEPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  thePhotoElectricEffect->AddEmModel(0, thePEPenelopeModel);

  //Compton scattering
  G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
  G4PenelopeComptonModel* theComptonPenelopeModel = new G4PenelopeComptonModel();
  theComptonScattering->SetEmModel(new G4KleinNishinaModel());
  theComptonPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  theComptonScattering->AddEmModel(0, theComptonPenelopeModel);

  //Gamma conversion
  G4GammaConversion* theGammaConversion = new G4GammaConversion();
  G4PenelopeGammaConversionModel* theGCPenelopeModel = new G4PenelopeGammaConversionModel();
  theGCPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  theGammaConversion->AddEmModel(0, theGCPenelopeModel);

  //Rayleigh scattering
  G4RayleighScattering* theRayleigh = new G4RayleighScattering();
  G4PenelopeRayleighModel* theRayleighPenelopeModel = new G4PenelopeRayleighModel();
  theRayleigh->SetEmModel(theRayleighPenelopeModel);

  ph->RegisterProcess(thePhotoElectricEffect, particle);
  ph->RegisterProcess(theComptonScattering, particle);
  ph->RegisterProcess(theGammaConversion, particle);
  ph->RegisterProcess(theRayleigh, particle);

  // e-
  particle = G4Electron::Electron();

  // multiple scattering
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
      
  //Ionisation
  G4eIonisation* eioni = new G4eIonisation();
  eioni->SetFluctModel(G4EmStandUtil::ModelOfFluctuations());
  G4PenelopeIonisationModel* theIoniPenelope = new G4PenelopeIonisationModel();
  theIoniPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);     
  eioni->AddEmModel(0, theIoniPenelope);
      
  //Bremsstrahlung
  G4eBremsstrahlung* brem = new G4eBremsstrahlung();
  G4PenelopeBremsstrahlungModel* theBremPenelope = new G4PenelopeBremsstrahlungModel();
  theBremPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  brem->SetEmModel(theBremPenelope);

  G4ePairProduction* ee = new G4ePairProduction();

  // register processes
  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(ss, particle);

  // e+
  particle = G4Positron::Positron();
    
  // multiple scattering
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
      
  //Ionisation
  eioni = new G4eIonisation();
  eioni->SetFluctModel(G4EmStandUtil::ModelOfFluctuations());
  theIoniPenelope = new G4PenelopeIonisationModel();
  theIoniPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  eioni->AddEmModel(0,theIoniPenelope);

  //Bremsstrahlung
  brem = new G4eBremsstrahlung();
  theBremPenelope = new G4PenelopeBremsstrahlungModel();
  theBremPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  brem->SetEmModel(theBremPenelope);
      
  //Annihilation
  G4eplusAnnihilation* anni = new G4eplusAnnihilation();
  G4PenelopeAnnihilationModel* theAnnPenelope = new G4PenelopeAnnihilationModel();
  theAnnPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
  anni->AddEmModel(0, theAnnPenelope);

  // register processes
  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(anni, particle);
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
