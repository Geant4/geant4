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
// $Id: TargetPhysicsList.cc,v 1.1 2003-10-08 12:32:11 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "TargetPhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

#include "G4LEProtonInelastic.hh"
#include "G4HadronKineticModel.hh"


TargetPhysicsList::TargetPhysicsList()
{;}

TargetPhysicsList::~TargetPhysicsList()
{;}

void TargetPhysicsList::ConstructParticle()
{
  G4Proton::ProtonDefinition();
}

void TargetPhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
  ConstructHad();
}

void TargetPhysicsList::ConstructHad()
{
  G4ProcessManager* pManager = 0;
  G4HadronKineticModel* theKineticModel = new G4HadronKineticModel();

  // proton
  //
  pManager = G4Proton::Proton()->GetProcessManager();
  G4LEProtonInelastic* theLEProtonModel = new G4LEProtonInelastic();
  theProtonInelastic.RegisterMe(theLEProtonModel);
  //  theProtonInelastic.RegisterMe(theKineticModel);
  pManager->AddDiscreteProcess(&theProtonInelastic);
}

void TargetPhysicsList::SetCuts()
{
  // suppress error messages even in case e/gamma/proton do not exist    
  G4int temp = GetVerboseLevel();
  SetVerboseLevel(0);                                    
  SetCutsWithDefault();   

  // Retrieve verbose level
  SetVerboseLevel(temp);  
}



