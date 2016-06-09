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
/// @file PhysicsList.cc
/// @brief Define physics list

#include "PhysicsList.hh"
#include "Particles.hh"
#include "PhysicsListEMstd.hh"
#include "G4SystemOfUnits.hh"

// --------------------------------------------------------------------------
PhysicsList::PhysicsList()
  :  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*mm;
  SetVerboseLevel(1);

  // particles
  RegisterPhysics(new Particles);

  // EM Physics
  RegisterPhysics(new PhysicsListEMstd);
}

// --------------------------------------------------------------------------
PhysicsList::~PhysicsList()
{
}

// --------------------------------------------------------------------------
void PhysicsList::SetCuts()
{
  SetCutsWithDefault();   
}
