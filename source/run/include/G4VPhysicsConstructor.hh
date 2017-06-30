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
// $Id: G4VPhysicsConstructor.hh 103953 2017-05-04 11:27:35Z gcosmo $
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
//    Add   G4PhysicsBuilderInterface	 21 Apr	 2017 by A.Dotti
// ------------------------------------------------------------
#ifndef G4VPhysicsConstructor_h
#define G4VPhysicsConstructor_h 1

#include "globals.hh"
#include "rundefs.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsListHelper.hh"
#include "G4VUPLSplitter.hh"
#include <vector>

class G4PhysicsBuilderInterface;

class G4VPCData
{
    //Encapsulate the fields of class G4VPhysicsConstructor
    //that are per-thread.
public:
    using PhysicsBuilders_V=std::vector<G4PhysicsBuilderInterface*>;
    void initialize();
    G4ParticleTable::G4PTblDicIterator* _aParticleIterator;
    PhysicsBuilders_V * _builders;
};

// The type G4VPCManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4VPCData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4VPhysicsConstructor instance,
// there is a corresponding G4VPCData instance. All G4VPCData instances
// are organized by the class G4VUPLManager as an array.
// The field "int g4vuplInstanceID" is added to the class G4VUserPhysicsList.
// The value of this field in each G4VUserPhysicsList instance is the
// subscript of the corresponding G44VUPLData instance.
// In order to use the class G44VUPLManager, we add a static member in the class
// G4VUserPhysicsList as follows: "static G4VUPLManager subInstanceManager".
// Both the master thread and worker threads change the length of the array
// for G44VUPLData instances mutually along with G4VUserPhysicsList
// instances are created. For each worker thread, it dynamically creates ions.
// Consider any thread A, if there is any other thread which creates an ion.
// This ion is shared by the thread A. So the thread A leaves an empty space
// in the array of G4PDefData instances for the ion.
//
// Important Note: you may wonder why we are introducing this mechanism
//                 since there is only one PL for each application.
//                 This is true, in the sense that only one PL is allowed
//                 to be associated to a G4RunManager, however user can instantiate
//                 as many PLs are needed and at run-time select one of the PLs to be used
//                 we thus need this mechanism to guarantee that the system works without
//                 problems in case of this (unusual) case. This may be reviewed in the future
typedef G4VUPLSplitter<G4VPCData> G4VPCManager;
typedef G4VPCManager G4VPhyscicsConstructorManager;

// This macros change the references to fields that are now encapsulated
// in the class G4VPCData.
//
// Note1:  the use of this-> this is needed to avoid compilation errors
// when using templated class with T=G4VUserPhysicsList. Don't know why.
// Note2: the name of the first #define is different, because otherwise
//        we need to change its use in all classes that inherits from
//        this base class (all examples). However one should note comment
//        on JIRA task: http://jira-geant4.kek.jp/browse/DEV-27

//#define aParticleIterator ((subInstanceManager.offset[g4vpcInstanceID])._aParticleIterator)

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
    G4int g4vpcInstanceID;
    G4RUN_DLL static G4VPCManager subInstanceManager;
    G4ParticleTable::G4PTblDicIterator* GetParticleIterator() const;
    using PhysicsBuilder_V=G4VPCData::PhysicsBuilders_V;
    //This returns a copy of the vector of pointers
    PhysicsBuilder_V GetBuilders() const;
    void AddBuilder(G4PhysicsBuilderInterface* bld);
public:
    inline G4int GetInstanceID() const;
    static const G4VPCManager& GetSubInstanceManager();

    //Method called by kernel to destroy thread-local
    //data, equivalent to destructor in sequential mode
    //Derived classes implementing this method, must also call
    //this base class method.
    virtual void TerminateWorker();
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
    return G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(process,particle);
  //return aPLHelper->RegisterProcess(process, particle);
}

inline
const G4VPCManager& G4VPhysicsConstructor::GetSubInstanceManager()
{
    return subInstanceManager;
}
#endif




