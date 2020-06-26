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

#include "G4EmLowEPPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// *** Processes and models

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LowEPComptonModel.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"

// e+-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"
#include "G4ePairProduction.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4Generator2BS.hh"
#include "G4BetheHeitler5DModel.hh"

// e+
#include "G4eplusAnnihilation.hh"

// hadrons
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4LindhardSorensenIonModel.hh"
#include "G4NuclearStopping.hh"

// msc models
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4LowEWentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4CoulombScattering.hh"

// interfaces
#include "G4LossTableManager.hh"
#include "G4EmBuilder.hh"
#include "G4EmParameters.hh"

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
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmLowEPPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLowEPPhysics::G4EmLowEPPhysics(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmLowEPPhysics"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  param->SetMinEnergy(100*CLHEP::eV);
  param->SetLowestElectronEnergy(100*CLHEP::eV);
  param->SetNumberOfBinsPerDecade(20);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetStepFunction(0.2, 100*CLHEP::um);
  param->SetStepFunctionMuHad(0.1, 50*CLHEP::um);
  param->SetStepFunctionLightIons(0.1, 20*CLHEP::um);
  param->SetStepFunctionIons(0.1, 1*CLHEP::um);
  param->SetUseMottCorrection(true);
  param->SetMscRangeFactor(0.04);   
  param->SetMuHadLateralDisplacement(true);
  param->SetFluo(true);
  param->SetUseICRU90Data(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmLowEPPhysics::~G4EmLowEPPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEPPhysics::ConstructParticle()
{
  // minimal set of particles for EM physics
  G4EmBuilder::ConstructMinimalEmSet();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmLowEPPhysics::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4EmBuilder::PrepareEMPhysics();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // common processes
  G4ePairProduction* ee = new G4ePairProduction();
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

  // nuclear stopping
  G4NuclearStopping* pnuc = new G4NuclearStopping();
  pnuc->SetMaxKinEnergy(CLHEP::MeV);

  // Add gamma EM Processes
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
  G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();
  pe->SetEmModel(theLivermorePEModel);

  // Compton scattering - Livermore model above 20 MeV, Monarsh's model below 
  G4ComptonScattering* cs = new G4ComptonScattering;
  cs->SetEmModel(new G4KleinNishinaModel());
  G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
  theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
  cs->AddEmModel(0, theLowEPComptonModel);

  // gamma conversion - 5D model below 80 GeV with Livermore x-sections
  G4GammaConversion* theGammaConversion = new G4GammaConversion();
  G4VEmModel* conv = new G4BetheHeitler5DModel();
  theGammaConversion->SetEmModel(conv);

  // default Rayleigh scattering is Livermore
  G4RayleighScattering* theRayleigh = new G4RayleighScattering();

  ph->RegisterProcess(pe, particle);
  ph->RegisterProcess(cs, particle);
  ph->RegisterProcess(theGammaConversion, particle);
  ph->RegisterProcess(theRayleigh, particle);

  // e-
  particle = G4Electron::Electron();

  // multiple scattering
  G4eMultipleScattering* msc = new G4eMultipleScattering();
  msc->SetEmModel(new G4LowEWentzelVIModel());

  // Ionisation - Livermore should be used only for low energies
  G4eIonisation* eioni = new G4eIonisation();
  G4LivermoreIonisationModel* theIoniLivermore = new G4LivermoreIonisationModel();
  theIoniLivermore->SetHighEnergyLimit(0.1*CLHEP::MeV); 
  eioni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );

  // Bremsstrahlung
  G4eBremsstrahlung* brem = new G4eBremsstrahlung();
  G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
  G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br1->SetHighEnergyLimit(GeV);
 
  // register processes
  ph->RegisterProcess(msc, particle);
  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);

  // e+
  particle = G4Positron::Positron();

  // multiple scattering
  msc = new G4eMultipleScattering();
  msc->SetEmModel(new G4LowEWentzelVIModel());

  // Standard ionisation 
  eioni = new G4eIonisation();

  // Bremsstrahlung
  brem = new G4eBremsstrahlung();
  br1 = new G4SeltzerBergerModel();
  br2 = new G4eBremsstrahlungRelModel();
  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br1->SetHighEnergyLimit(CLHEP::GeV);

  // register processes
  ph->RegisterProcess(msc, particle);
  ph->RegisterProcess(eioni, particle);
  ph->RegisterProcess(brem, particle);
  ph->RegisterProcess(ee, particle);
  ph->RegisterProcess(new G4eplusAnnihilation(), particle);

  // generic ion
  particle = G4GenericIon::GenericIon();
       
  G4ionIonisation* ionIoni = new G4ionIonisation();
  G4VEmModel* mod = new G4LindhardSorensenIonModel();
  ionIoni->SetEmModel(mod);

  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);
  ph->RegisterProcess(pnuc, particle);

  // muons, hadrons ions
  G4EmBuilder::ConstructCharged(hmsc, pnuc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
