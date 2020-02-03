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
// ====================================================================
//   ExN01PhysicsList.cc
//
//                                         2005 Q
// ====================================================================
#include "ExN01PhysicsList.hh"
#include "G4ParticleTypes.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////////////
ExN01PhysicsList::ExN01PhysicsList()
////////////////////////////////////
{
}

/////////////////////////////////////
ExN01PhysicsList::~ExN01PhysicsList()
/////////////////////////////////////
{
}

//////////////////////////////////////////
void ExN01PhysicsList::ConstructParticle()
//////////////////////////////////////////
{
  G4Geantino::GeantinoDefinition();
}

/////////////////////////////////////////
void ExN01PhysicsList::ConstructProcess()
/////////////////////////////////////////
{
  AddTransportation();
}

////////////////////////////////
void ExN01PhysicsList::SetCuts()
////////////////////////////////
{
  G4int temp = GetVerboseLevel();
  SetVerboseLevel(0);

  SetCutsWithDefault();   

  SetVerboseLevel(temp);  
}

