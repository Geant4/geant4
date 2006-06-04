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
// $Id: ExN01PhysicsList.cc,v 1.3 2006-06-04 21:36:35 kmura Exp $
// $Name: not supported by cvs2svn $
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

