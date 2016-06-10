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
// $Id: G4GenericPhysicsList.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: Witek Pokorski
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef TG4GenericPhysicsList_h
#define TG4GenericPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

#include "G4PhysicsConstructorRegistry.hh"

#include "G4GenericMessenger.hh"

template<class T>
class TG4GenericPhysicsList: public T
{
public:
  TG4GenericPhysicsList(G4int ver = 1);
  TG4GenericPhysicsList(std::vector<G4String>* physConstr, G4int ver = 1);
  virtual ~TG4GenericPhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

  void RegisterPhysicsConstructor(G4String& physconstr) {this->RegisterPhysics(G4PhysicsConstructorRegistry::Instance()->GetPhysicsConstructor(physconstr));}

  G4GenericMessenger messenger;
  void DeclareProperties();

};
#include "G4GenericPhysicsList.icc"
typedef TG4GenericPhysicsList<G4VModularPhysicsList> G4GenericPhysicsList;

#endif



