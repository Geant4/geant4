//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VModularPhysicsList.hh,v 1.4 2002-05-29 12:11:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//   In addition the user must implement the four virtual methods
//   in his own concrete class derived from this class. 
//        G4VModularPhysicsList::SetCuts()
//   to set cut values in range to all particles
//   (rebuild of physics table will be invoked).

// ------------------------------------------------------------ 
// History
// - first version                   12 Nov 2000 by H.Kurashige 
// ------------------------------------------------------------
#ifndef G4VModularPhysicsList_h
#define G4VModularPhysicsList_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"

#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"

class G4VModularPhysicsList: public virtual G4VUserPhysicsList
{
  public: 
    G4VModularPhysicsList();
    virtual ~G4VModularPhysicsList();

  public:  // with description
    //  "SetCuts" method sets a cut value for all particle types 
    //   in the particle table
    virtual void SetCuts() = 0; 

  protected: // with description
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();
  
  public: // with description
    void RegisterPhysics(G4VPhysicsConstructor* );

    const G4VPhysicsConstructor* GetPhysics(G4int index) const;
    const G4VPhysicsConstructor* GetPhysics(const G4String& name) const;

  /////////////////////////////////////

  protected: // with description
   // vector of pointers to G4VPhysicsConstructor
   typedef G4std::vector<G4VPhysicsConstructor*> G4PhysConstVector;
   G4PhysConstVector* physicsVector;
};
   

inline 
 void G4VModularPhysicsList::RegisterPhysics(G4VPhysicsConstructor* fPhysics)
{
  physicsVector->push_back(fPhysics);
}    

inline  
 const G4VPhysicsConstructor* G4VModularPhysicsList::GetPhysics(G4int idx) const
{
  G4int i;
  G4PhysConstVector::iterator itr= physicsVector->begin();
  for (i=0; i<idx && itr!= physicsVector->end() ; ++i) ++itr;
  if (itr!= physicsVector->end()) return (*itr);
  else return 0;
}

inline  
 const G4VPhysicsConstructor* G4VModularPhysicsList::GetPhysics(const G4String& name) const
{
  G4PhysConstVector::iterator itr;
  for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
    if ( name == (*itr)->GetPhysicsName()) break;
  }
  if (itr!= physicsVector->end()) return (*itr);
  else return 0;
}

   
    

#endif
