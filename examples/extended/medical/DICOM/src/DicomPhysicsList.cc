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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//
// History: 30.11.07  P.Arce default cut changed to 1 mm
//*******************************************************

#include "DicomPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

DicomPhysicsList::DicomPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1.*mm;
  cutForGamma     = 1.*mm;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  SetVerboseLevel(0);
}

DicomPhysicsList::~DicomPhysicsList()
{}

void DicomPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructBosons();
  ConstructLeptons();
}

void DicomPhysicsList::ConstructBosons()
{
  // gamma
  G4Gamma::GammaDefinition();

}

void DicomPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void DicomPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

#include "G4MultipleScattering.hh"
// gamma
#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
// e-
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
// e+
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

void DicomPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      //processes
      lowePhot = new  G4LowEnergyPhotoElectric("LowEnPhotoElec");
      loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
      loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");

      if (particleName == "gamma")
        {
	  //gamma
	  pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
	  pmanager->AddDiscreteProcess(lowePhot);
	  pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
	  pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);

        }
      else if (particleName == "e-")
        {
	  //electron
	  pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
	  pmanager->AddProcess(loweIon,     -1, 2,2);
	  pmanager->AddProcess(loweBrem,    -1,-1,3);

        }
      else if (particleName == "e+")
        {
	  //positron
	  pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
	  pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
	  pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);

        }
    }
}
#include "G4Decay.hh"
void DicomPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (theDecayProcess->IsApplicable(*particle))
        {
	  pmanager ->AddProcess(theDecayProcess);
	  // set ordering for PostStepDoIt and AtRestDoIt
	  pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
	  pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
        }
    }
}

void DicomPhysicsList::SetCuts()
{
  if (verboseLevel >0)
    {
      G4cout << "DicomPhysicsList::SetCuts:";
      G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");


  if (verboseLevel>0)
    DumpCutValuesTable();
}

void DicomPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

void DicomPhysicsList::SetElectronCut(G4double val)
{
  //  ResetCuts();
  cutForElectron = val;
}

void DicomPhysicsList::SetPositronCut(G4double val)
{
  //  ResetCuts();
  cutForPositron = val;
}

