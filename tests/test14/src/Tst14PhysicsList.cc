//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: Tst14PhysicsList.cc,v 1.17 2003-02-23 14:36:31 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Redisegned for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst14PhysicsList.hh"
#include "Tst14PhysicsListMessenger.hh"
#include "Tst14DetectorConstruction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"

#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4LowEnergyPolarizedCompton.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

// e+
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

Tst14PhysicsList::Tst14PhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false)
{
  defaultCutValue = 0.1 * mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;

  SetVerboseLevel(1);

  // UI messenger
  physicsListMessenger = new Tst14PhysicsListMessenger(this);
 
  // Particles
  RegisterPhysics( new Tst14Particles("particles") );

  // General chunk of PhysicsList (transportation,...)
  RegisterPhysics( new Tst14GeneralProcesses("general") );
}


Tst14PhysicsList::~Tst14PhysicsList()
{
  delete physicsListMessenger;
}

void Tst14PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

void Tst14PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      // gamma    
      if (comptonLowEPolarised)
	{
	  G4std::cout << "Loading Polarised Compton" << G4std::endl;
	  pmanager->AddDiscreteProcess(new G4LowEnergyPolarizedCompton);
	}
      else
	{
	  G4std::cout << "Loading LowEnergy Compton" << G4std::endl;
	  pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
	}
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);

      LePeprocess = new G4LowEnergyPhotoElectric();
      pmanager->AddDiscreteProcess(LePeprocess);

      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);

      LeIoprocess = new G4LowEnergyIonisation();
      pmanager->AddProcess(LeIoprocess, -1,  2, 2);

      LeBrprocess = new G4LowEnergyBremsstrahlung();
      pmanager->AddProcess(LeBrprocess, -1, -1, 3);

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      
    } 
  }
}

void Tst14PhysicsList::SetGELowLimit(G4double cut)
{
  if (verboseLevel > 0)
    {
      G4cout << "Tst14PhysicsList - Gamma and electron cut in energy = " 
	     << cut * MeV << " MeV" << G4endl;
    }  
  G4Gamma::SetEnergyRange(cut,1e5);
  G4Electron::SetEnergyRange(cut,1e5);
  G4Positron::SetEnergyRange(cut,1e5);
}

void Tst14PhysicsList::SetGammaLowLimit(G4double cut)
{
  if (verboseLevel > 0)
    {
      G4cout << "Tst14PhysicsList - Gamma cut in energy = " 
	     << cut * MeV << " MeV" << G4endl;
    }  
  G4Gamma::SetEnergyRange(cut,1e5);
}

void Tst14PhysicsList::SetElectronLowLimit(G4double cut)
{
  if (verboseLevel > 0)
    {
      G4cout << "Tst14PhysicsList - Electron cut in energy = " 
	     << cut * MeV << " MeV" << G4endl;
    }  
  G4Electron::SetEnergyRange(cut,1e5);
}

void Tst14PhysicsList::SetGammaCut(G4double value)
{
  ResetCuts();
  cutForGamma = value;
}


void Tst14PhysicsList::SetElectronCut(G4double value)
{
  ResetCuts();
  cutForElectron = value;
}

void Tst14PhysicsList::SetCuts()
{
  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
}


void Tst14PhysicsList::SetLowEnSecPhotCut(G4double cut)
{
  // This m.f. is pertinent to LowEnergy EPDL/EEDL processes only

  // Get the ProcessManager for photons and the list of photon processes
  G4ProcessManager* photonManager = G4Gamma::GammaDefinition()->GetProcessManager();
  G4ProcessVector* photonProcesses = photonManager->GetProcessList();
  G4int nPhotonProcesses = photonProcesses->size();

  // Loop over photon processes until one retrieves LowEnergyPhotoElectric
  for (G4int iPhoton=0; i<nPhotonProcesses; i++)
    {
      G4VProcess* process = photonProcesses[iPhoton];
      G4String& name = process->GetProcessName();
      G4String nameLowE("LowEnPhotoElec");
      if (name == nameLowE)
	{
	  // The only way to get access to the cut stting is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyPhotoElectric* lowEProcess = dynamic_cast<G4LowEnergyPhotoElectric*> process;
	  if (lowEProcess != 0) 
	    {
	     lowEProcess->SetCutForLowEnSecPhotons(cut);
	     G4cout << "Low energy secondary photons cut is now set to: "
		    << cut * MeV
		    << " (MeV) for LowEnergyPhotoElectric"
		    << G4endl;
	    }
	}
    }
  
  // Get the ProcessManager for electrons and the list of electron processes
  G4ProcessManager* electronManager = G4Electron::ElectronDefinition()->GetProcessManager();
  G4ProcessVector* electronProcesses = electronManager->GetProcessList();
  G4int nElectronProcesses = electronProcesses->size();

  // Loop over electron processes until one retrieves LowEnergyIonisation or LowEnergyBremsstrahlung
  for (G4int iElectron=0; i<nElectronProcesses; i++)
    {
      G4VProcess* process = electronProcesses[iElectron];
      G4String& name = process->GetProcessName();

      G4String nameIoni("LowEnergyIoni");
      if (name == nameIoni)
	{
	  // The only way to get access to the cut setting is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyIonisation* lowEProcess = dynamic_cast<G4LowEnergyIonisation*> process;
	  if (lowEProcess != 0) 
	    {
	      lowEProcess->SetCutForLowEnSecPhotons(cut);
	      G4cout << "Low energy secondary photons cut is now set to: "
		     << cut * MeV
		     << " (MeV) for LowEnergyIonisation"
		     << G4endl;
	    }
	}

      G4String nameBrems("LowEnBrem");
      if (name == nameBrems)
	{
	  // The only way to get access to the cut setting is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyIonisation* lowEProcess = dynamic_cast<G4LowEnergyBremsstrahlung*> process;
	  if (lowEProcess != 0) 
	    {
	      lowEProcess->SetCutForLowEnSecPhotons(cut);
	      G4cout << "Low energy secondary photons cut is now set to: "
		     << cut * MeV
		     << " (MeV) for LowEnergyBremsstrahlung"
		     << G4endl;
	    }
	}
    }
}

void Tst14PhysicsList::SetLowEnSecElecCut(G4double cut)
{  
  // This m.f. is pertinent to LowEnergy EPDL/EEDL processes only

  // Get the ProcessManager for photons and the list of photon processes
  G4ProcessManager* photonManager = G4Gamma::GammaDefinition()->GetProcessManager();
  G4ProcessVector* photonProcesses = photonManager->GetProcessList();
  G4int nPhotonProcesses = photonProcesses->size();

  // Loop over photon processes until one retrieves LowEnergyPhotoElectric
  for (G4int iPhoton=0; i<nPhotonProcesses; i++)
    {
      G4VProcess* process = photonProcesses[iPhoton];
      G4String& name = process->GetProcessName();
      G4String nameLowE("LowEnPhotoElec");
      if (name == nameLowE)
	{
	  // The only way to get access to the cut stting is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyPhotoElectric* lowEProcess = dynamic_cast<G4LowEnergyPhotoElectric*> process;
	  if (lowEProcess != 0) 
	    {
	     lowEProcess->SetCutForLowEnSecElectrons(cut);
	     G4cout << "Low energy secondary electrons cut is now set to: "
		    << cut * MeV
		    << " (MeV) for LowEnergyPhotoElectric"
		    << G4endl;
	    }
	}
    }
  
  // Get the ProcessManager for electrons and the list of electron processes
  G4ProcessManager* electronManager = G4Electron::ElectronDefinition()->GetProcessManager();
  G4ProcessVector* electronProcesses = electronManager->GetProcessList();
  G4int nElectronProcesses = electronProcesses->size();

  // Loop over electron processes until one retrieves LowEnergyIonisation
  for (G4int iElectron=0; i<nElectronProcesses; i++)
    {
      G4VProcess* process = electronProcesses[iElectron];
      G4String& name = process->GetProcessName();
      G4String nameLowE("LowEnergyIoni");
      if (name == nameLowE)
	{
	  // The only way to get access to the cut setting is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyIonisation* lowEProcess = dynamic_cast<G4LowEnergyIonisation*> process;
	  if (lowEProcess != 0) 
	    {
	      lowEProcess->SetCutForLowEnSecElectrons(cut);
	      G4cout << "Low energy secondary electrons cut is now set to: "
		     << cut * MeV
		     << " (MeV) for LowEnergyIonisation"
		     << G4endl;
	    }
	}
    }
}

void Tst14PhysicsList::ActivateAuger(G4bool value)
{  
  // Get the ProcessManager for photons and the list of photon processes
  G4ProcessManager* photonManager = G4Gamma::GammaDefinition()->GetProcessManager();
  G4ProcessVector* photonProcesses = photonManager->GetProcessList();
  G4int nPhotonProcesses = photonProcesses->size();

  // Loop over photon processes until one retrieves LowEnergyPhotoElectric
  for (G4int iPhoton=0; i<nPhotonProcesses; i++)
    {
      G4VProcess* process = photonProcesses[iPhoton];
      G4String& name = process->GetProcessName();
      G4String nameLowE("LowEnPhotoElec");
      if (name == nameLowE)
	{
	  // The only way to get access to the Auger activation is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyPhotoElectric* lowEProcess = dynamic_cast<G4LowEnergyPhotoElectric*> process;
	  if (lowEProcess != 0) 
	    {
	     lowEProcess->ActivateAuger(value);
	      G4cout << "Auger electron production flag is " << value
		     << " for LowEnergyPhotoElectric process" << G4endl;
	    }
	}
    }
  
  // Get the ProcessManager for electrons and the list of electron processes
  G4ProcessManager* electronManager = G4Electron::ElectronDefinition()->GetProcessManager();
  G4ProcessVector* electronProcesses = electronManager->GetProcessList();
  G4int nElectronProcesses = electronProcesses->size();

  // Loop over electron processes until one retrieves LowEnergyIonisation
  for (G4int iElectron=0; i<nElectronProcesses; i++)
    {
      G4VProcess* process = electronProcesses[iElectron];
      G4String& name = process->GetProcessName();
      G4String nameLowE("LowEnergyIoni");
      if (name == nameLowE)
	{
	  // The only way to get access to the Auger activation is through a dynamic_cast
	  // (it is ugly!)
	  G4LowEnergyIonisation* lowEProcess = dynamic_cast<G4LowEnergyIonisation*> process;
	  if (lowEProcess != 0) 
	    {
	      lowEProcess->ActivateAuger(value);
	      G4cout << "Auger electron production flag is " << value
		     << " for LowEnergyIonisation process" << G4endl;
	    }
	}
    }
}








