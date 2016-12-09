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
#ifndef G4PhysicsConstructorFactory_h
#define G4PhysicsConstructorFactory_h 1


#include "globals.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4VPhysicsConstructor.hh"

class G4VBasePhysConstrFactory
{

public:

  virtual G4VPhysicsConstructor* Instantiate() = 0;

};


template <typename T> class G4PhysicsConstructorFactory : public G4VBasePhysConstrFactory
{
public:
  
  G4PhysicsConstructorFactory(const G4String& name)
  {
    G4PhysicsConstructorRegistry::Instance()->AddFactory(name, this);
  }
  
  virtual G4VPhysicsConstructor* Instantiate() 
  {
    return new T();
  }
};


#define G4_DECLARE_PHYSCONSTR_FACTORY(physics_constructor) \
  const G4PhysicsConstructorFactory<physics_constructor>& physics_constructor##Factory = G4PhysicsConstructorFactory<physics_constructor>(#physics_constructor)

// support for phys constructors defined within a namespace
// a bit tricky due to cpp macro expansion and the use of "::"
// use  G4_DECLARE_PHYSCONSTR_FACTORY_NS( myns::MyProc, myns, MyProc );
#define G4_DECLARE_PHYSCONSTR_FACTORY_NS( physics_constructor, nsname, pcbase )  \
  namespace nsname { \
    const G4PhysicsConstructorFactory<physics_constructor>& pcbase##Factory = G4PhysicsConstructorFactory<physics_constructor>(#physics_constructor); \
  } \
  typedef int xyzzy__LINE__
  // eat trailing semicolon using silly typedef

// REFERENCE (rather than DECLARE) when the physics constructor if it is part
// of a static library.  No need to include the header (DECLARE needs this
// to build the code), we just need to make a reference in order to pull
// the compilation unit static variable from the library and cause it
// to be initialized (and thus self-register)

#define G4_REFERENCE_PHYSCONSTR_FACTORY(physics_constructor) \
  class physics_constructor; \
  extern const G4PhysicsConstructorFactory<physics_constructor>& physics_constructor##Factory; \
  const G4PhysicsConstructorFactory<physics_constructor>& physics_constructor##FactoryRef = physics_constructor##Factory

#define G4_REFERENCE_PHYSCONSTR_FACTORY_NS( physics_constructor, nsname, pcbase ) \
  namespace nsname { \
    class pcbase; \
    extern const G4PhysicsConstructorFactory<physics_constructor>& pcbase##Factory; \
    const G4PhysicsConstructorFactory<physics_constructor>& pcbase##FactoryRef = pcbase##Factory; \
  } \
  typedef int xyzzy__LINE__
  // eat trailing semicolon using silly typedef

#endif
