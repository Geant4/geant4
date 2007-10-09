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
// -------------------------------------------------------------------
// $Id: MicrodosimetryPhysicsList.cc,v 1.1 2007-10-09 08:00:29 sincerti Exp $
// -------------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4StepLimiter.hh"
#include "G4BaryonConstructor.hh"	              
#include "G4IonConstructor.hh"	 
#include "G4MesonConstructor.hh"	 

#include "MicrodosimetryPhysicsList.hh"
//DNA
#include "G4DNAGenericIonsManager.hh"
// END DNA

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetryPhysicsList::MicrodosimetryPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 0.01*micrometer;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetryPhysicsList::~MicrodosimetryPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructBaryons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPhysicsList::ConstructBosons()
{ 
  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPhysicsList::ConstructBaryons()
{
  //  baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();
  //DNA

  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");

//END DNA

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructHad();
  ConstructGeneral();
}

#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyRayleigh.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

#include "G4hLowEnergyIonisation.hh"
#include "G4hMultipleScattering.hh"

// DNA
#include "G4VLowEnergyTestableDiscreteProcess.hh"

#include "G4DNAScreenedRutherford.hh"
#include "G4DNAExcitationEmfietzoglou.hh"

#include "G4DNAIonizationBorn.hh"
#include "G4DNAExcitationMillerGreen.hh"
#include "G4DNAExcitationBorn.hh"
#include "G4DNAIonizationRudd.hh"
#include "G4DNAChargeDecreaseDingfelder.hh"
#include "G4DNAChargeIncreaseDingfelder.hh"
//END DNA


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrodosimetryPhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {


      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);     

      G4LowEnergyPhotoElectric * LePeprocess = new G4LowEnergyPhotoElectric();
      LePeprocess->ActivateAuger(true);
      LePeprocess->SetCutForLowEnSecPhotons(0.250 * keV);
      LePeprocess->SetCutForLowEnSecElectrons(0.250 * keV);
      pmanager->AddDiscreteProcess(LePeprocess);
      
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());
      
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh());
      
      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 3);

      
    } else if (particleName == "e-") {

// DNA

       pmanager->AddDiscreteProcess(new G4DNAScreenedRutherford());
       pmanager->AddDiscreteProcess(new G4DNAExcitationEmfietzoglou());
       pmanager->AddDiscreteProcess(new G4DNAIonizationBorn());

// END DNA

        } else if (particleName == "e+") {

      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      
      pmanager->AddProcess(new G4eIonisation,      -1, 2,2);
      
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      
      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 3);

    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

    } else if ( particleName == "alpha" || 
                particleName == "alpha+" || 
		particleName == "helium" ||
		particleName == "proton" ||
		particleName == "hydrogen"  
		 ) {
     
// DNA
       pmanager->AddDiscreteProcess(new G4DNAIonizationRudd());
       pmanager->AddDiscreteProcess(new G4DNAExcitationMillerGreen());
       pmanager->AddDiscreteProcess(new G4DNAIonizationBorn());
       pmanager->AddDiscreteProcess(new G4DNAExcitationBorn());
       pmanager->AddDiscreteProcess(new G4DNAChargeDecreaseDingfelder());
       pmanager->AddDiscreteProcess(new G4DNAChargeIncreaseDingfelder());
       pmanager->AddProcess(new G4StepLimiter(),-1,-1,6);       
      }

// END DNA

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4LEAlphaInelastic.hh"

void MicrodosimetryPhysicsList::ConstructHad()
{

  G4HadronElasticProcess * theElasticProcess = new G4HadronElasticProcess;
  theElasticProcess->RegisterMe( new G4LElastic() );

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pManager = particle->GetProcessManager();

    if (particle->GetParticleName() == "alpha") 
       { 

 	  // INELASTIC SCATTERING
	  // Binary Cascade
  	  G4BinaryLightIonReaction* theBC = new G4BinaryLightIonReaction();
  	  theBC -> SetMinEnergy(80.*MeV);
  	  theBC -> SetMaxEnergy(40.*GeV);
  
  	  // TRIPATHI CROSS SECTION
  	  // Implementation of formulas in analogy to NASA technical paper 3621 by 
  	  // Tripathi, et al. Cross-sections for ion ion scattering
  	  G4TripathiCrossSection* TripathiCrossSection = new G4TripathiCrossSection;
  
  	  // IONS SHEN CROSS SECTION
  	  // Implementation of formulas 
  	  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  	  // Total Reaction Cross Section for Heavy-Ion Collisions
  	  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;
  
  	  // Final state production model for Alpha inelastic scattering below 20 GeV
  	  G4LEAlphaInelastic* theAIModel = new G4LEAlphaInelastic;
  	  theAIModel -> SetMaxEnergy(100.*MeV);

	  G4AlphaInelasticProcess * theIPalpha = new G4AlphaInelasticProcess;		  
	  theIPalpha->AddDataSet(TripathiCrossSection);
	  theIPalpha->AddDataSet(aShen);

	  // Register the Alpha Inelastic and Binary Cascade Model
	  theIPalpha->RegisterMe(theAIModel);
	  theIPalpha->RegisterMe(theBC);
	  
	  // Activate the alpha inelastic scattering using the alpha inelastic and binary cascade model
	  pManager -> AddDiscreteProcess(theIPalpha);
	  
	  // Activate the Hadron Elastic Process
	  pManager -> AddDiscreteProcess(theElasticProcess); 
            
       }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MicrodosimetryPhysicsList::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "MicrodosimetryPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  if (verboseLevel>0) DumpCutValuesTable();

}

