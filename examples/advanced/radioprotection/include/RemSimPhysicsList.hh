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
// $Id$
//
// Author: Susanna Guatelli, susanna@uow.edu.au

#ifndef REMSIMPHYSICSLIST_HH
#define REMSIMPHYSICSLIST_HH 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "G4EmConfigurator.hh"

class G4VPhysicsConstructor;

class RemSimPhysicsListMessenger;

class RemSimPhysicsList: public G4VModularPhysicsList {
public:
  
  RemSimPhysicsList();

  virtual ~RemSimPhysicsList();

  void ConstructParticle();

  void ConstructProcess();

  //register the threshold of production of secondaries
  void SetCuts();
  
  // Register PhysicsList chunks
  void AddPhysicsList(const G4String& name);

private:
  G4EmConfigurator em_config;
  G4bool helIsRegistered;// hadronic elastic scattering flag
  G4bool bicIsRegistered;// binary cascade inleastic scattering flag
  G4bool bicIonIsRegistered; // binary ion cascade inelastic scattering flag
  G4bool radioactiveDecayIsRegistered;// radioactive decay module flag

  G4VPhysicsConstructor*               emPhysicsList;
  G4VPhysicsConstructor*               decPhysicsList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;

  RemSimPhysicsListMessenger* messenger;
};
#endif







