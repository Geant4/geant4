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

// ====================================================================
//
//   ExN04PrimaryGeneratorAction.hh
//   $Id: ExN04PrimaryGeneratorAction.hh,v 1.1 2002-04-29 20:43:45 asaim Exp $
//
// ====================================================================
#ifndef EXN04_PRIMARY_GENERATOR_ACTION_H
#define EXN04_PRIMARY_GENERATOR_ACTION_H

#include "g4std/map"
#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class G4VPrimaryGenerator;
class ExN04PrimaryGeneratorMessenger;

class ExN04PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
  G4VPrimaryGenerator* particleGun;
  G4VPrimaryGenerator* hepmcAscii;
  G4VPrimaryGenerator* pythiaGen;

  G4VPrimaryGenerator* currentGenerator;
  G4String currentGeneratorName;
  std::map<G4String, G4VPrimaryGenerator*> gentypeMap;

  ExN04PrimaryGeneratorMessenger* messenger;

public:
  ExN04PrimaryGeneratorAction();
  ~ExN04PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* anEvent);

  void SetGenerator(G4VPrimaryGenerator* gen);
  void SetGenerator(G4String genname);

  G4VPrimaryGenerator* GetGenerator() const;
  G4String GetGeneratorName() const;
};

// ====================================================================
// inline functions
// ====================================================================

inline void ExN04PrimaryGeneratorAction::SetGenerator
(G4VPrimaryGenerator* gen) 
{ 
  currentGenerator= gen; 
}

inline void ExN04PrimaryGeneratorAction::SetGenerator(G4String genname) 
{
  map<G4String, G4VPrimaryGenerator*>::iterator pos= gentypeMap.find(genname);
  if(pos != gentypeMap.end()) {
    currentGenerator= pos->second;
    currentGeneratorName= genname;
  }  
}

inline G4VPrimaryGenerator* ExN04PrimaryGeneratorAction::GetGenerator() const 
{ 
  return currentGenerator; 
}

inline G4String ExN04PrimaryGeneratorAction::GetGeneratorName() const 
{ 
  return currentGeneratorName; 
}

#endif
