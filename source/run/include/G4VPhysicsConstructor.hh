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
// $Id$
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
//      Only one physics constructor can be registered to
//      Modular Physics List for each "physics_type".
//      Physics constructors with same "physics_type" can be
//      replaced by using the method of 
//      G4VModularPhysicsList::ReplacePhysics()
//      
//
// ------------------------------------------- 
//   History
//    first version                      12 Nov. 2000 by H.Kurashige 
//    Add   physicsType                  14 Mar. 2011 by H.Kurashige
//    Add   RegisterProcess               1 May  2011 by H.Kurashige
// ------------------------------------------------------------
#ifndef G4VPhysicsConstructor_h
#define G4VPhysicsConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsListHelper.hh"

class G4VPhysicsConstructor
{
  public:  // with description

    G4VPhysicsConstructor(const G4String& ="");
    G4VPhysicsConstructor(const G4String& name, G4int physics_type);
    virtual ~G4VPhysicsConstructor();

    virtual void ConstructParticle()=0;
      // This method will be invoked in the Construct() method. 
      // each particle type will be instantiated
 
    virtual void ConstructProcess()=0;
      // This method will be invoked in the Construct() method.
      // each physics process will be instantiated and
      // registered to the process manager of each particle type 

    inline void  SetPhysicsName(const G4String& ="");
    inline const G4String& GetPhysicsName() const;

    inline void  SetPhysicsType(G4int);
    inline G4int GetPhysicsType() const;

    inline void  SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
      // set/get controle flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More
      // verbose level is set equal to physics list when registered 

  protected: 

    inline G4bool RegisterProcess(G4VProcess*            process,
                                  G4ParticleDefinition*  particle);
      // Register a process to the particle type 
      // according to the ordering parameter table
      // 'true' is returned if the process is registerd successfully

  protected:
    G4int    verboseLevel;
    G4String namePhysics;
    G4int    typePhysics;

    G4ParticleTable* theParticleTable;
    G4ParticleTable::G4PTblDicIterator* theParticleIterator;
      // the particle table has the complete List of existing particle types

    G4PhysicsListHelper* thePLHelper;
};

// Inlined methods

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

inline void G4VPhysicsConstructor::SetPhysicsType(G4int val)
{
  if (val>0) typePhysics = val;
}

inline G4int G4VPhysicsConstructor::GetPhysicsType() const
{
  return  typePhysics;
}

inline 
 G4bool G4VPhysicsConstructor::RegisterProcess(G4VProcess*            process,
					       G4ParticleDefinition*  particle)
{
  return thePLHelper->RegisterProcess(process, particle);
}
#endif




