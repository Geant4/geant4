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
// $Id: Tst22PhysicsList.cc,v 1.1 2001-11-15 15:10:13 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst22PhysicsList.hh"

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

#include "Tst22GeneralPhysics.hh"
#include "Tst22EMPhysics.hh"
#include "Tst22MuonPhysics.hh"
#include "Tst22HadronPhysics.hh"
#include "Tst22IonPhysics.hh"

Tst22PhysicsList::Tst22PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new Tst22GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new Tst22EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new Tst22MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new Tst22HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new Tst22IonPhysics("ion"));

}

Tst22PhysicsList::~Tst22PhysicsList()
{
}

void Tst22PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}



