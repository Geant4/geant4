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
// $Id: Tst23PhysicsList.cc,v 1.3 2004-03-05 15:25:45 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst23PhysicsList.hh"

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
#include <iomanip>   

#include "Tst23GeneralPhysics.hh"
#include "Tst23EMPhysics.hh"
#include "Tst23MuonPhysics.hh"
#include "Tst23HadronPhysics.hh"
#include "Tst23IonPhysics.hh"

Tst23PhysicsList::Tst23PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new Tst23GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new Tst23EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new Tst23MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new Tst23HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new Tst23IonPhysics("ion"));

}

Tst23PhysicsList::~Tst23PhysicsList()
{
}

void Tst23PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}



