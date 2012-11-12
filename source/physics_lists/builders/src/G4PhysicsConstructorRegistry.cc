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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4PhysicsConstructorRegistry
//
// Author  W. Pokorski  21.09.2012
//
// Modifications:
//

#include "G4ios.hh"

#include "G4PhysicsConstructorRegistry.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4PhysicsConstructorFactory.hh"

G4PhysicsConstructorRegistry* G4PhysicsConstructorRegistry::theInstance = 0;

G4PhysicsConstructorRegistry* G4PhysicsConstructorRegistry::Instance()
{
  if(0 == theInstance) {
    static G4PhysicsConstructorRegistry manager;
    theInstance = &manager;
  }
  return theInstance;
}

G4PhysicsConstructorRegistry::G4PhysicsConstructorRegistry()
{}

G4PhysicsConstructorRegistry::~G4PhysicsConstructorRegistry()
{
  Clean();
}

void G4PhysicsConstructorRegistry::Clean()
{
  size_t n = physConstr.size(); 
  if(n > 0) {
    for (size_t i=0; i<n; ++i) {
      if(physConstr[i]) {
	G4VPhysicsConstructor* p = physConstr[i];
	physConstr[i] = 0;
	delete p;
      }
    }
    physConstr.clear();
  }
}

void G4PhysicsConstructorRegistry::Register(G4VPhysicsConstructor* p)
{
  if(!p) return;
  size_t n = physConstr.size(); 
  if(n > 0) {
    for (size_t i=0; i<n; ++i) {
      if(physConstr[i] == p) { return; }
    }
  }
  physConstr.push_back(p);
}

void G4PhysicsConstructorRegistry::DeRegister(G4VPhysicsConstructor* p)
{
  if(!p) return;
  size_t n = physConstr.size(); 
  if(n > 0) {
    for (size_t i=0; i<n; ++i) {
      if(physConstr[i] == p) {
        physConstr[i] = 0;
	return;
      }
    }
  }
}

void G4PhysicsConstructorRegistry::AddFactory(G4String name, G4VBasePhysConstrFactory* factory)
{
  factories[name] = factory;
}

G4VPhysicsConstructor* G4PhysicsConstructorRegistry::GetPhysicsConstructor(const G4String& name)
{
  // check if factory exists...
  //
  if (factories.find(name)!=factories.end())
    {
        // we could store the list of called factories in some vector and
        // before returning we can could first check if this physics constructor was already instantiated
        // if yes, we can throw an exception saying that this physics can been already registered
        
      return factories[name]->Instantiate();
    }
  else
    {
      G4ExceptionDescription ED;
      ED << "The factory for the physics constructor ["<< name << "] does not exist!" << G4endl;
      G4Exception("G4PhysicsConstructorRegistry::GetPhysicsConstructor", "PhysicsList001", FatalException, ED);
      return 0;
    }
}
