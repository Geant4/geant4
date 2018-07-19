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
// $Id: G4VModularPhysicsList.hh 103803 2017-04-27 14:03:05Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// Class Description:
//   This class is a subclass of G4VUserPhysicsList.     
//   The user should register his/her physics constructors 
//   by using 
//         G4VModularPhysicsList::RegsiterPhysics() 
//   to construt particles and processes.
//
//   Only one physics constructor can be registered for each "physics_type".
//   Physics constructors with same "physics_type" can be replaced by
//   G4VModularPhysicsList::ReplacePhysics() method
//
// ------------------------------------------------------------ 
// History
// - first version                   12 Nov 2000 by H.Kurashige 
// - Add  ReplacePhysics             14 Mar 2011 by H.Kurashige
// - Add Worker cleanup              21 Apr 2017 by A.Dotti
//
// ------------------------------------------------------------
#ifndef G4VModularPhysicsList_h
#define G4VModularPhysicsList_h 1

#include <vector>

#include "globals.hh"
#include "rundefs.hh"
#include "G4ios.hh"

#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUPLSplitter.hh"

class G4VMPLData {
    //Encapsulate the fields of class G4VModularPhysicsList
    //that are per-thread.
public:
    void initialize();
    typedef std::vector<G4VPhysicsConstructor*> G4PhysConstVectorData;
    //TODO: understand this
    //See: https://jira-geant4.kek.jp/browse/DEV-284
    G4PhysConstVectorData* physicsVector;
};

// The type G4VMPLManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4VMPLData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4VUserPhysicsList instance,
// there is a corresponding G4VMPLData instance. All G4VMPLData instances
// are organized by the class G4VMPLManager as an array.
// The field "int G4VMPLInstanceID" is added to the class G4VUserPhysicsList.
// The value of this field in each G4VUserPhysicsList instance is the
// subscript of the corresponding G44VUPLData instance.
// In order to use the class G44VUPLManager, we add a static member in the class
// G4VUserPhysicsList as follows: "static G4VMPLManager subInstanceManager".
// Both the master thread and worker threads change the length of the array
// for G44VUPLData instances mutually along with G4VUserPhysicsList
// instances are created.
typedef G4VUPLSplitter<G4VMPLData> G4VMPLManager;
typedef G4VMPLManager G4VModularPhysicsListSubInstanceManager;

class G4VModularPhysicsList: public virtual G4VUserPhysicsList
{
  public: 
    G4VModularPhysicsList();
    virtual ~G4VModularPhysicsList();
  
 protected:
    // hide copy constructor and assignment operator
    G4VModularPhysicsList(const G4VModularPhysicsList&);
    G4VModularPhysicsList & operator=(const G4VModularPhysicsList&);

  public:  // with description
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle() override;
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess() override;
  
  public: // with description
    // Register Physics Constructor 
    void RegisterPhysics(G4VPhysicsConstructor* );
    
    const G4VPhysicsConstructor* GetPhysics(G4int index) const;
    const G4VPhysicsConstructor* GetPhysics(const G4String& name) const;
    const G4VPhysicsConstructor* GetPhysicsWithType(G4int physics_type) const;

    // Replace Physics Constructor 
    //  The existing physics constructor with same physics_type as one of
    //  the given physics constructor is replaced
    //  (existing physics will be deleted)
    //  If a corresponding physics constructor is NOT found, 
    //  the given physics constructor is just added         
    void ReplacePhysics(G4VPhysicsConstructor* );

   // Remove Physics Constructor from the list
    void RemovePhysics(G4VPhysicsConstructor* );
    void RemovePhysics(G4int type);
    void RemovePhysics(const G4String& name);
    
  /////////////////////////////////////
  public: // with description
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // set/get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More
   // given verbose level is set to all physics constructors 

  protected: // with description
   G4int verboseLevel;
    typedef G4VMPLData::G4PhysConstVectorData G4PhysConstVector;
    G4int g4vmplInstanceID;
    G4RUN_DLL static G4VMPLManager G4VMPLsubInstanceManager;
  public:
    inline G4int GetInstanceID() const;
    static const G4VMPLManager& GetSubInstanceManager();
    virtual void TerminateWorker() override;
};
   
inline  
 G4int G4VModularPhysicsList::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline
G4int G4VModularPhysicsList::GetInstanceID() const
{
    return g4vmplInstanceID;
}

inline
const G4VMPLManager& G4VModularPhysicsList::GetSubInstanceManager()
{
    return G4VMPLsubInstanceManager;
}
#endif
