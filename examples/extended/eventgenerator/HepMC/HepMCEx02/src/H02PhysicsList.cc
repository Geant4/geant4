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
// $Id: H02PhysicsList.cc,v 1.1 2002-05-28 14:15:48 murakami Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "H02PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "H02GeneralPhysics.hh"
#include "H02EMPhysics.hh"
#include "H02MuonPhysics.hh"
#include "H02HadronPhysics.hh"
#include "H02IonPhysics.hh"

H02PhysicsList::H02PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new H02GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new H02EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new H02MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new H02HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new H02IonPhysics("ion"));


}

H02PhysicsList::~H02PhysicsList()
{
}

void H02PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}



