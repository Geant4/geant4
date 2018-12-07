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
//

#include "globals.hh"

#include "eRositaPhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4hImpactIonisation.hh"

#include "G4ProductionCutsTable.hh"


eRositaPhysicsList::eRositaPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 0.001*mm;
   SetVerboseLevel(1);

   std::cout << "==============================================================================="
	     << std::endl
	     << "Geant4 eRosita example - based on a simplified version of eROSITA simulation"
	     << std::endl
	     << "Further details can be found in:"
	     << std::endl
	     << "M.G. Pia et al., 'PIXE Simulation With Geant4', "
	     << "IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, 2009"
	     << std::endl
	     << "N. Meidinger et al., 'Development of the focal plane PNCCD camera system for the X-ray space telescope eROSITA', " 
	     << std::endl
	     <<"NIM A 624, 321-329, 2010"
	     << std::endl
	     << "==============================================================================="
	     << std::endl;

   std::cout<< std::endl;
   
   std::cout << "==============================================================================="
	     << std::endl
	     << " The use of G4LowEnergyIonisation, G4LowEnergyBremsstrahlung, "
	     << std::endl
	     << "G4LowEnergyPhotoElectric, G4LowEnergyCompton, G4LowEnergyGammaConversion"
	     << std::endl
	     << "in this example is intentional. These classes will be replaced by other classes"
	     << std::endl
	     << "appropriate to the problem domain in a forthcoming Geant4 version"
	     << std::endl
	     << "==============================================================================="
	     << std::endl;
}


eRositaPhysicsList::~eRositaPhysicsList()
{}


void eRositaPhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
}


void eRositaPhysicsList::ConstructBosons()
{
  // pseudo-particles
  //G4Geantino::GeantinoDefinition();
  //G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}


void eRositaPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  //G4MuonPlus::MuonPlusDefinition();
  //G4MuonMinus::MuonMinusDefinition();
  // nu_e
  //G4NeutrinoE::NeutrinoEDefinition();
  //G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  //G4NeutrinoMu::NeutrinoMuDefinition();
  //G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}


void eRositaPhysicsList::ConstructMesons()
{
  //  mesons
  //    light mesons
  //G4PionPlus::PionPlusDefinition();
  //G4PionMinus::PionMinusDefinition();
  //G4PionZero::PionZeroDefinition();
  //G4Eta::EtaDefinition();
  //G4EtaPrime::EtaPrimeDefinition();
  //G4KaonPlus::KaonPlusDefinition();
  //G4KaonMinus::KaonMinusDefinition();
  //G4KaonZero::KaonZeroDefinition();
  //G4AntiKaonZero::AntiKaonZeroDefinition();
  //G4KaonZeroLong::KaonZeroLongDefinition();
  //G4KaonZeroShort::KaonZeroShortDefinition();
}


void eRositaPhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  //G4Neutron::NeutronDefinition();
  //G4AntiNeutron::AntiNeutronDefinition();
}


void eRositaPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
  //AddStepMax();
}



void eRositaPhysicsList::ConstructEM()
{
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* processManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      // photon   

      G4PhotoElectricEffect* photoelectric = new G4PhotoElectricEffect;
      //photoelectric->ActivateAuger(true);
      //photoelectric->SetCutForLowEnSecPhotons(0.250 * keV);
      //photoelectric->SetCutForLowEnSecElectrons(0.250 * keV);
      G4ComptonScattering* compton = new G4ComptonScattering;
      G4GammaConversion* gammaConversion = new G4GammaConversion;
      G4RayleighScattering* rayleigh = new G4RayleighScattering;

      processManager -> AddDiscreteProcess(photoelectric);
      processManager -> AddDiscreteProcess(compton);
      processManager -> AddDiscreteProcess(gammaConversion);
      processManager -> AddDiscreteProcess(rayleigh);
    
    } else if (particleName == "e-") {

      // electron

      G4eMultipleScattering* eMultipleScattering = new G4eMultipleScattering();
      G4eIonisation* eIonisation = new G4eIonisation();
      G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();

      processManager -> AddProcess(eMultipleScattering, -1, 1, 1);
      processManager -> AddProcess(eIonisation, -1, 2, 2);
      processManager -> AddProcess(eBremsstrahlung, -1, -1, 3);   

    } else if (particleName == "e+") {
      // positron
      processManager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      processManager->AddProcess(new G4eIonisation,         -1, 2, 2);
      processManager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
      processManager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);

      //} else if( particleName == "mu+" || 
      //        particleName == "mu-"    ) {
      //muon  
      //processManager->AddProcess(new G4MuMultipleScattering, -1, 1, 1);
      //processManager->AddProcess(new G4MuIonisation,         -1, 2, 2);
      //processManager->AddProcess(new G4MuBremsstrahlung,     -1, 3, 3);
      //processManager->AddProcess(new G4MuPairProduction,     -1, 4, 4);       
             
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {
      //proton  
      /*
      G4hImpactIonisation* hIonisation = new G4hImpactIonisation();
      hIonisation->SetPixeCrossSectionK("ecpssr");
      hIonisation->SetPixeCrossSectionL("ecpssr");
      hIonisation->SetPixeCrossSectionM("ecpssr");
      hIonisation->SetPixeProjectileMinEnergy(1.* keV);
      hIonisation->SetPixeProjectileMaxEnergy(200. * MeV);
      hIonisation->SetCutForSecondaryPhotons(250. * eV);
      hIonisation->SetCutForAugerElectrons(250. * eV);
      */
      G4hIonisation* hIonisation = new G4hIonisation();

      G4hMultipleScattering* hMultipleScattering = new G4hMultipleScattering();

      processManager -> AddProcess(hMultipleScattering, -1, 1, 1);   
      processManager -> AddProcess(hIonisation, -1, 2, 2);
     
    } else if( particleName == "alpha" || 
	       particleName == "He3" ||
	       particleName == "pi-" ||
               particleName == "pi+" ||
	       particleName == "GenericIon" ) {

      // pions, alpha, ions (should never occur in the current example) 
      processManager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      processManager->AddProcess(new G4ionIonisation,       -1, 2, 2);
                      
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      processManager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      processManager->AddProcess(new G4hIonisation,         -1, 2, 2);
    }
  }
}

#include "G4Decay.hh"

void eRositaPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* processManager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      processManager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      processManager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      processManager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}
  

/*
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void eRositaPhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  G4StepLimiter* stepLimiter = new G4StepLimiter();
  ////G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();
  
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* processManager = particle->GetProcessManager();

      if (particle->GetPDGCharge() != 0.0)
        {
	  processManager ->AddDiscreteProcess(stepLimiter);
	  ////processManager ->AddDiscreteProcess(userCuts);
        }
  }
}
*/

void eRositaPhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 
  //
  SetCutsWithDefault();

  // Set the secondary production cut lower than 990. eV
  // Very important for processes at low energies
 
  G4double lowLimit = 250. * eV;
  G4double highLimit = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
     
  if (verboseLevel>0) DumpCutValuesTable();
}


