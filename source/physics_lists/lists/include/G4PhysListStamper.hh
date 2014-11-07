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
#ifndef G4PhysListStamper_h
#define G4PhysListStamper_h 1


#include "globals.hh"
#include "G4PhysListRegistry.hh"
#include "G4VModularPhysicsList.hh"

class G4VBasePhysListStamper
{

public:

  virtual G4VModularPhysicsList* Instantiate(G4int /* verbose */) = 0;

};


template <typename T> class G4PhysListStamper : public G4VBasePhysListStamper
{
public:
  
  G4PhysListStamper(const G4String& name)
  {
    G4PhysListRegistry::Instance()->AddFactory(name, this);
  }
  
  virtual G4VModularPhysicsList* Instantiate(G4int verbose) 
  {
    return new T(verbose);
  }
};

#define G4_DECLARE_PHYSLIST_FACTORY(physics_list) \
  const G4PhysListStamper<physics_list>& physics_list##Factory = G4PhysListStamper<physics_list>(#physics_list)

// support for physics list defined within a namespace
// a bit tricky because cpp macro expansion doesn't like "::"
// ala  G4_DECLARE_PHYSLIST_FACTORY_NS(myns::MyPL,myns,MyPL)   // without trailing ";"
#define G4_DECLARE_PHYSLIST_FACTORY_NS( physics_list , nsname , plbase )    \
  namespace nsname {                                                    \
    const G4PhysListStamper<physics_list>& plbase##Factory = G4PhysListStamper<physics_list>(#physics_list); \
  }

#endif
