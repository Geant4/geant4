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
// $Id: 
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4PhysicsConstructorRegistry
//
// Author  W. Pokorski  21.09.2012
//
// Modifications:
// 
//
// Class Description
// This is a singleton keeping pointers to all physics constructors
// Class Description - End

#ifndef G4PhysicsConstructorRegistry_h
#define G4PhysicsConstructorRegistry_h 1

#include <vector>
#include <map>
#include "globals.hh"


class G4VPhysicsConstructor;
class G4VBasePhysConstrFactory;

class G4PhysicsConstructorRegistry
{
public:

  static G4PhysicsConstructorRegistry* Instance();
  // access 
  
  ~G4PhysicsConstructorRegistry();
  
  void Register(G4VPhysicsConstructor*);
  //register new physics constructors

  void DeRegister(G4VPhysicsConstructor*);
  //deregister physics constructors

  void Clean();
  //clean the store
  
  void AddFactory(G4String, G4VBasePhysConstrFactory*);

  G4VPhysicsConstructor* GetPhysicsConstructor(const G4String& name);

  G4bool IsKnownPhysicsConstructor(const G4String& name);

  std::vector<G4String> AvailablePhysicsConstructors() const;

  void PrintAvailablePhysicsConstructors() const;
    
private:

  G4PhysicsConstructorRegistry();

  static G4ThreadLocal G4PhysicsConstructorRegistry* theInstance;
  
  std::vector <G4VPhysicsConstructor*> physConstr;

  std::map <G4String, G4VBasePhysConstrFactory*> factories;

};

#endif
