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
// $Id: DicomPhysicsList.cc 75679 2013-11-05 09:08:41Z gcosmo $
//
/// \file medical/DICOM/src/DicomPhysicsList.cc
/// \brief Implementation of the DicomPhysicsList class
//

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4StepLimiter.hh"
#include "G4BaryonConstructor.hh"                      
#include "G4IonConstructor.hh"         
#include "G4MesonConstructor.hh"         
#include "G4SystemOfUnits.hh"

#include "DicomPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DicomPhysicsList::DicomPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 0.01*micrometer;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DicomPhysicsList::~DicomPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomPhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructBaryons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomPhysicsList::ConstructBosons()
{ 
  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomPhysicsList::ConstructBaryons()
{
  //  baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructHad();
  ConstructGeneral();
}

// *** Processes and models

// gamma

#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

// e-

#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+

#include "G4eplusAnnihilation.hh"

// mu

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadrons

#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

//

#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomPhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {


      G4PhotoElectricEffect* thePhotoElectricEffect =
                                                   new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel
                                         = new G4LivermorePhotoElectricModel();
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel =
                                                 new G4LivermoreComptonModel();
      theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
      pmanager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel
                                       = new G4LivermoreGammaConversionModel();
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      pmanager->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel
                                              = new G4LivermoreRayleighModel();
      theRayleigh->AddEmModel(0, theRayleighModel);
      pmanager->AddDiscreteProcess(theRayleigh);
      
      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 5);
      
    } else if (particleName == "e-") {
     
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      
      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->AddEmModel(0, new G4LivermoreIonisationModel(),
                           new G4UniversalFluctuation() );
      eIoni->SetStepFunction(0.2, 100*um); //     
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      
      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->AddEmModel(0, new G4LivermoreBremsstrahlungModel());
      pmanager->AddProcess(eBrem,                   -1,-3, 3);

      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 4);
      
    } else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
    
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);
      
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);
      
      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 5);

    } else if (particleName == "GenericIon") {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 20*um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);

      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 3);

    } else if (particleName == "alpha" ||
               particleName == "He3" ) {

      // Identical to G4EmStandardPhysics_option3
      
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 20*um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);

      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 3);
    }

   //

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"

void DicomPhysicsList::ConstructHad()
{

  G4HadronElasticProcess * theElasticProcess = new G4HadronElasticProcess;
  theElasticProcess->RegisterMe( new G4HadronElastic() );

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pManager = particle->GetProcessManager();

    if (particle->GetParticleName() == "alpha") { 

      // INELASTIC SCATTERING
      // Binary Cascade
      G4BinaryLightIonReaction* theBC = new G4BinaryLightIonReaction();
      // DHW - change lower limit from 80 MeV to zero because
      //       G4LEDeuteronInelastic is deprecated
      theBC->SetMinEnergy(0.*MeV);
      theBC->SetMaxEnergy(40.*GeV);
  
      // TRIPATHI CROSS SECTION
      // Implementation of formulas in analogy to NASA technical paper 3621 by 
      // Tripathi, et al. Cross-sections for ion ion scattering
      G4TripathiCrossSection* TripathiCrossSection = new G4TripathiCrossSection;
  
      // IONS SHEN CROSS SECTION
      // Implementation of formulas 
      // Shen et al. Nuc. Phys. A 491 130 (1989) 
      // Total Reaction Cross Section for Heavy-Ion Collisions
      G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;
  
      G4AlphaInelasticProcess* theIPalpha = new G4AlphaInelasticProcess;                  
      theIPalpha->AddDataSet(TripathiCrossSection);
      theIPalpha->AddDataSet(aShen);

      // Register the Binary Cascade Model
      theIPalpha->RegisterMe(theBC);
          
      // Activate the alpha inelastic scattering using the binary cascade model
      pManager -> AddDiscreteProcess(theIPalpha);
          
      // Activate the Hadron Elastic Process
      pManager -> AddDiscreteProcess(theElasticProcess);           
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomPhysicsList::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DicomPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "DicomPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(fCutForGamma, "gamma");
  SetCutValue(fCutForElectron, "e-");
  SetCutValue(fCutForPositron, "e+");
  
  if (verboseLevel>0) DumpCutValuesTable();

}

