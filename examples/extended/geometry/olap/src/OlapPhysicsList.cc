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
//
// $Id: OlapPhysicsList.cc,v 1.1 2002-06-04 07:40:22 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapPhysicsList
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapPhysicsList.hh"

#include "G4ios.hh"              
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"

OlapPhysicsList::OlapPhysicsList():  G4VUserPhysicsList()
{
  currentDefaultCut = defaultCutValue = 2.0*mm;
  cutForGamma       = defaultCutValue;
  cutForElectron    = defaultCutValue;
  cutForProton      = defaultCutValue;

 SetVerboseLevel(1);
}


OlapPhysicsList::~OlapPhysicsList()
{}


void OlapPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
}


void OlapPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
}


void OlapPhysicsList::ConstructProcess()
{
  AddTransportation();
  
  theParticleIterator->reset();
  //theParticleTable->DumpTable();
  
  return;
  
  // Olap-Process for stopping at user-defined surface ...    
/*  will come in the next release
  while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();      
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddProcess(new OlapProc(),-1,-1,1);
  }
*/
}


#include "G4Decay.hh"


void OlapPhysicsList::SetCuts()
{
  // don't need cuts because transportation is the only process ..
  return; 
}
