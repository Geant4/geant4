// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPhysicsConstructor.hh,v 1.1 2000-11-14 23:53:13 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an virtual class for constructing 
//      particles and processes. This class objects will be
//      registered to G4VPhysicsList.
//
//      User must implement following four virtual methods
//      in his own concrete class derived from this class. 
//
//      all necessary particle type will be instantiated
//      virtual void ConstructParticle();
// 
//      all physics processes will be instantiated and
//      registered to the process manager of each particle type 
//      virtual void ConstructProcess();
//
//
// ------------------------------------------- 
//	History
//        first version                   12 Nov. 2000 by H.Kurashige 
// ------------------------------------------------------------
#ifndef G4VPhysicsConstructor_h
#define G4VPhysicsConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
class G4VPhysicsConstructor
{
  public: 
    G4VPhysicsConstructor(const G4String& ="");
    virtual ~G4VPhysicsConstructor(){};

  public: // with description
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle()=0;
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess()=0;

  public: // with description
   void  SetPhysicsName(const G4String& ="");
   const G4String& GetPhysicsName() const;

   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // set/get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More
   // verbose level is set equal to physics list when registered 

 protected:
   G4String namePhysics;
   G4int verboseLevel;
 
 protected:
  // the particle table has the complete List of existing particle types
  G4ParticleTable* theParticleTable;
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;

};


inline 
 G4VPhysicsConstructor::G4VPhysicsConstructor(const G4String& name)
   :verboseLevel(0),namePhysics(name)
{
  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
}

inline void G4VPhysicsConstructor::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline  G4int G4VPhysicsConstructor::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline void G4VPhysicsConstructor::SetPhysicsName(const G4String& name)
{
  namePhysics = name;
}

inline const G4String&  G4VPhysicsConstructor::GetPhysicsName() const
{
  return  namePhysics;
}

#endif




