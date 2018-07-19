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
// $Id: G4NuclearPolarizationStore.hh 96490 2016-04-19 06:57:04Z gcosmo $
//
// Created:
// 03 July 2017 V.Ivanchenko
//
// Class Description
// This is a thread local singleton class to store all G4NuclearPolarization
// objects
//
// Class Description - End

#ifndef G4NuclearPolarizationStore_h
#define G4NuclearPolarizationStore_h 1

#include "globals.hh"
#include "G4NuclearPolarization.hh"
#include "G4ThreadLocalSingleton.hh"

const G4int maxNumStates = 10;

class G4NuclearPolarizationStore
{
friend class G4ThreadLocalSingleton<G4NuclearPolarizationStore>;

public:

  static G4NuclearPolarizationStore* GetInstance();
  // access 
  
  ~G4NuclearPolarizationStore();

  void Register(G4NuclearPolarization* ptr);
  // register new G4NuclearPolarization object

  G4NuclearPolarization* FindOrBuild(G4int Z, G4int A, G4double Eexc);
  // find G4NuclearPolarization object or build new

  void RemoveMe(G4NuclearPolarization* ptr);
  // deregister and delete G4NuclearPolarization object

private:

  G4NuclearPolarizationStore();

  static G4ThreadLocal G4NuclearPolarizationStore* instance;

  G4NuclearPolarization* nuclist[maxNumStates];
  G4int oldIdx;
};

#endif
