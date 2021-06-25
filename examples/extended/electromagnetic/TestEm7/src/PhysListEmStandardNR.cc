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
/// \file electromagnetic/TestEm7/src/PhysListEmStandardNR.cc
/// \brief Implementation of the PhysListEmStandardNR class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandardNR.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4UrbanMscModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4Generator2BS.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

#include "G4ScreenedNuclearRecoil.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardNR::PhysListEmStandardNR(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetStepFunction(0.2, 100*um);
  param->SetStepFunctionMuHad(0.2, 50*um);          
  param->SetStepFunctionLightIons(0.1, 10*um);
  param->SetStepFunctionIons(0.1, 1*um);
  param->SetFluo(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardNR::~PhysListEmStandardNR()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandardNR::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // muon & hadron bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
  G4MuPairProduction* mup = new G4MuPairProduction();

  G4ScreenedNuclearRecoil* nucr = new G4ScreenedNuclearRecoil();
  G4double energyLimit = 100.*MeV;
  nucr->SetMaxEnergyForScattering(energyLimit);
  G4eCoulombScatteringModel* csm = new G4eCoulombScatteringModel();
  csm->SetActivationLowEnergyLimit(energyLimit);

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

       // Compton scattering
      G4ComptonScattering* cs = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel());
      ph->RegisterProcess(cs, particle);

      // Photoelectric
      G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
      G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();
      theLivermorePEModel->SetHighEnergyLimit(10*GeV);
      pe->SetEmModel(theLivermorePEModel);
      ph->RegisterProcess(pe, particle);

      // Gamma conversion
      G4GammaConversion* gc = new G4GammaConversion();
      G4VEmModel* thePenelopeGCModel = new G4PenelopeGammaConversionModel();
      thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
      gc->SetEmModel(thePenelopeGCModel);
      ph->RegisterProcess(gc, particle);

      // Rayleigh scattering
      ph->RegisterProcess(new G4RayleighScattering(), particle);
      
    } else if (particleName == "e-") {

      // ionisation
      G4eIonisation* eIoni = new G4eIonisation();

      // bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();

      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(eBrem, particle);
            
    } else if (particleName == "e+") {
      // ionisation
      G4eIonisation* eIoni = new G4eIonisation();

      // bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();

      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(eBrem, particle);

      // annihilation at rest and in flight
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
            
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      G4MuIonisation* muIoni = new G4MuIonisation();

      ph->RegisterProcess(muIoni, particle);
      ph->RegisterProcess(mub, particle);
      ph->RegisterProcess(mup, particle);
      ph->RegisterProcess(new G4CoulombScattering(), particle);
             
    } else if (particleName == "alpha" || particleName == "He3") {

      G4hMultipleScattering* msc = new G4hMultipleScattering();
      G4UrbanMscModel* model = new G4UrbanMscModel();
      model->SetActivationLowEnergyLimit(energyLimit);
      msc->SetEmModel(model);
      ph->RegisterProcess(msc, particle);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ph->RegisterProcess(ionIoni, particle);

      ph->RegisterProcess(nucr, particle);

    } else if (particleName == "GenericIon" ) { 

      G4hMultipleScattering* msc = new G4hMultipleScattering();
      G4UrbanMscModel* model = new G4UrbanMscModel();
      model->SetActivationLowEnergyLimit(energyLimit);
      msc->SetEmModel(model);
      ph->RegisterProcess(msc, particle);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ph->RegisterProcess(ionIoni, particle);

      ph->RegisterProcess(nucr, particle);
     
    } else if (particleName == "proton" ||
               particleName == "deuteron" ||
               particleName == "triton") { 

      G4hMultipleScattering* msc = new G4hMultipleScattering();
      G4UrbanMscModel* model = new G4UrbanMscModel();
      model->SetActivationLowEnergyLimit(energyLimit);
      msc->SetEmModel(model);
      ph->RegisterProcess(msc, particle);

      G4hIonisation* hIoni = new G4hIonisation();
      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(nucr, particle);

    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4hIonisation(), particle);
    }
  }
  
  // Deexcitation
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

