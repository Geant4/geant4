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
// $Id: OlapPhysicsList.cc,v 1.3 2010-08-16 08:23:55 kurasige Exp $
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
