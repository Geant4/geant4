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
// $Id: G4VModularPhysicsList.cc,v 1.1 2002-05-29 03:48:00 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an derived class of G4VUserPhysicsList.     
//       User should regsiter his/her physics constructors 
//      by using 
//         G4VModularPhysicsList::RegsiterPhysics() 
//      to construt particles and processes.
//       In addition User must implement following four virtual methods
//      in his own concrete class derived from this class. 
//        G4VModularPhysicsList::SetCuts()
//           set cut values in range to all particles
//           (and rebuilding physics table will be invoked )
//
// ------------------------------------------- 
//	History
//        first version                   12 Nov. 2000 by H.Kurashige 
// ------------------------------------------------------------
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

#include "G4VModularPhysicsList.hh"


G4VModularPhysicsList::G4VModularPhysicsList()
                  : G4VUserPhysicsList()
{
   physicsVector = new G4PhysConstVector();
}

G4VModularPhysicsList::~G4VModularPhysicsList()
{
  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    delete (*itr);
  }
  physicsVector->clear();
}

void G4VModularPhysicsList::ConstructParticle()
{
  // create particles
  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    (*itr)->ConstructParticle();;
  }
}


void G4VModularPhysicsList::ConstructProcess()
{
  AddTransportation();

  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    (*itr)->ConstructProcess();
  }
}

