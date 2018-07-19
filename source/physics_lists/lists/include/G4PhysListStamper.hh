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
// ala  G4_DECLARE_PHYSLIST_FACTORY_NS(myns::MyPL,myns,MyPL);
#define G4_DECLARE_PHYSLIST_FACTORY_NS( physics_list , nsname , plbase )    \
  namespace nsname {                                                    \
    const G4PhysListStamper<physics_list>& plbase##Factory = G4PhysListStamper<physics_list>(#physics_list); \
  } \
  typedef int xyzzy__LINE__
  // eat trailing semicolon using silly typedef

// REFERENCE (rather than DECLARE) when the physics list if it is part
// of a static library.  No need to include the header (DECLARE needs this
// to build the code), we just need to make a reference in order to pull
// the compilation unit static variable from the library and cause it
// to be initialized (and thus self-register)

// this is _very_ much complicated because of the templating of these lists

#define G4_REFERENCE_PHYSLIST_FACTORY(physics_list)	\
  class G4VModularPhysicsList;			        \
  template <class T> class T##physics_list;	\
  typedef T##physics_list<G4VModularPhysicsList> physics_list;	       \
  extern const G4PhysListStamper<physics_list>& physics_list##Factory; \
  const G4PhysListStamper<physics_list>& physics_list##FactoryRef = physics_list##Factory

#define G4_REFERENCE_PHYSLIST_FACTORY_NS(physics_list, nsname, plbase )	\
  class G4VModularPhysicsList;    \
  namespace nsname { \
    template <class T> class T##plbase;	\
    typedef T##plbase<G4VModularPhysicsList> plbase; \
    extern const G4PhysListStamper<plbase>& plbase##Factory; \
    const G4PhysListStamper<plbase>& plbase##FactoryRef = plbase##Factory; \
  } \
  typedef int xyzzy__LINE__
  // eat trailing semicolon using silly typedef

#endif

