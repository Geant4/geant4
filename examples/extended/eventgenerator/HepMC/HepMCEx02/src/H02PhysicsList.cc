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
// $Id: H02PhysicsList.cc,v 1.4 2006/06/29 17:11:21 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
#include <iomanip>   

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



