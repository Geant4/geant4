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
// $Id: Sc01PhysicsList.cc,v 1.2 2004-11-10 07:43:14 grichine Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              10  Jan. 1998 by H.Kurashige
//        add decay at rest          26  Feb. 1998 by H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "Sc01PhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Geantino.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>                


Sc01PhysicsList::Sc01PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Sc01PhysicsList::~Sc01PhysicsList()
{
}

void Sc01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();

}

void Sc01PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}



void Sc01PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}
#include "G4OpBoundaryProcess.hh"

void Sc01PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "opticalphoton") 
    {
    // opticalphoton
    // Construct processes for opticalphoton
    //  pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());
 
    } 
  }
}

#include "G4Decay.hh"
void Sc01PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (theDecayProcess->IsApplicable(*particle)) 
    { 
      pmanager->AddProcess(theDecayProcess, INT_MAX, -1, INT_MAX); 
    }
  }
}

void Sc01PhysicsList::SetCuts()
{
  if (verboseLevel >0)
  {
    G4cout << "Sc01PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }
  // set verbose level 0 to surpress messages
  G4int temp = GetVerboseLevel();
  SetVerboseLevel(0);  

  SetCutsWithDefault();   

  // retrieve verbose level
  SetVerboseLevel(temp);
}


