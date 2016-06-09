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
// $Id: PhysicsList.cc,v 1.1 2007-11-16 14:16:46 kmura Exp $
// $Name: not supported by cvs2svn $
//
// ====================================================================
//   PhysicsList.cc
//
//                                         2007 Q
// ====================================================================
#include "PhysicsList.hh"
#include "Particles.hh"
#include "PhysicsListEMstd.hh"

// ====================================================================
//
// class description
//
// ====================================================================

///////////////////////////
PhysicsList::PhysicsList()
  :  G4VModularPhysicsList()
////////////////////////////
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*mm;
  SetVerboseLevel(1);

  // particles
  RegisterPhysics(new Particles);

  // EM Physics
  RegisterPhysics(new PhysicsListEMstd);
}


///////////////////////////
PhysicsList::~PhysicsList()
///////////////////////////
{
}


///////////////////////////
void PhysicsList::SetCuts()
///////////////////////////
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

