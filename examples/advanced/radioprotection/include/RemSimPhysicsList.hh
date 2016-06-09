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
// $Id: RemSimPhysicsList.hh,v 1.5 2004/05/22 12:57:05 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it

#ifndef REMSIMPHYSICSLIST_HH
#define REMSIMPHYSICSLIST_HH 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class RemSimPhysicsListMessenger;

class RemSimPhysicsList: public G4VModularPhysicsList {
public:
  
  RemSimPhysicsList();

  virtual ~RemSimPhysicsList();

  //register the threshold of production of secondaries
  virtual void SetCuts();
  
  // Register PhysicsList chunks
  void AddPhysicsList(const G4String& name);

private:
  G4bool electronIsRegistered;
  G4bool positronIsRegistered;
  G4bool photonIsRegistered;
  G4bool ionIsRegistered;
  G4bool hadronicIsRegistered; 
  G4bool decayIsRegistered;
  G4bool muonIsRegistered;
  RemSimPhysicsListMessenger* messenger;
};
#endif







