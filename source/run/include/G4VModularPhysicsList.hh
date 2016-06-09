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
//
// ------------------------------------------------------------
#ifndef G4VModularPhysicsList_h
#define G4VModularPhysicsList_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"

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
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();
  
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
    //  If any corresponding physics constructor is found, 
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
   // vector of pointers to G4VPhysicsConstructor
   typedef std::vector<G4VPhysicsConstructor*> G4PhysConstVector;
   G4PhysConstVector* physicsVector;
   G4int verboseLevel;
};
   


inline  
 G4int G4VModularPhysicsList::GetVerboseLevel() const
{
  return  verboseLevel;
}   
    

#endif
