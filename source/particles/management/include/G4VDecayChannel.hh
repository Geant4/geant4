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
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      30 May  1997 H.Kurashige
//      23 Mar. 2000 H.Weber      : add GetAngularMomentum()
// ------------------------------------------------------------
#ifndef G4VDecayChannel_h
#define G4VDecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include <cmath>

class    G4ParticleDefinition;
class    G4DecayProducts;
class    G4ParticleTable;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The class DecayChannelPrivateSubclass is introduced to
//encapsulate the fields associated to the class G4VDecayChannel that
//may not be read-only.
#ifndef DECAYCHANNELPRIVATESUBCLASS_HH
#define DECAYCHANNELPRIVATESUBCLASS_HH

class DecayChannelPrivateSubclass
{
public:
  G4ParticleDefinition*  parent;
  G4ParticleDefinition** daughters;
  G4double               parent_mass;
  G4double*              daughters_mass;
  void initialize() {
    parent = 0;
    daughters = 0;
    parent_mass = 0.0;
    daughters_mass = 0;
  };
};
#endif

//01.25.2009 Xin Dong: Phase II change for Geant4 multithreading.
//The class G4DecayChannelSubInstanceManager is introduced to 
//encapsulate the methods used by both the master thread and 
//worker threads to allocate memory space for the fields encapsulated
//by the class DecayChannelPrivateSubclass. When each thread
//initializes the value for these fields, it refers to them using a macro
//definition defined below. For every G4DecayChannel instance, there is
//a corresponding DecayChannelPrivateSubclass instance. All
//DecayChannelPrivateSubclass instances are organized by the
//class G4DecayChannelSubInstanceManager as an array. The field "  
//int g4decayChannelSubInstanceID" is added to the class G4DecayChannel.
//The value of this field in each G4DecayChannel instance is the subscript
//of the corresponding DecayChannelPrivateSubclass instance. In order
//to use the class G4DecayChannelSubInstanceManager, we add a static member in
//the class G4DecayChannel as follows: "  
//static G4DecayChannelSubInstanceManager g4decayChannelSubInstanceManager".
//For the master thread, the array for DecayChannelPrivateSubclass 
//instances grows dynamically along with G4DecayChannel instances are
//created. For each worker thread, it copies the array of 
//DecayChannelPrivateSubclass instances from the master thread.
//In addition, it invokes a method similiar to the constructor explicitly
//to achieve the partial effect for each instance in the array.

#ifndef G4DECAYCHANNELSUBINSTANCEMANAGER_HH
#define G4DECAYCHANNELSUBINSTANCEMANAGER_HH

#include "G4MTTransitoryParticle.hh"
typedef G4MTPrivateParticleCounter<DecayChannelPrivateSubclass> G4DecayChannelSubInstanceManager;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//These macros changes the references to fields that are now encapsulated
//in the class DecayChannelPrivateSubclass.
#define parentG4MTThreadPrivate ((G4VDecayChannel::g4decayChannelSubInstanceManager.offset[g4decayChannelInstanceID]).parent)
#define daughtersG4MTThreadPrivate ((G4VDecayChannel::g4decayChannelSubInstanceManager.offset[g4decayChannelInstanceID]).daughters)
#define parent_mass ((G4VDecayChannel::g4decayChannelSubInstanceManager.offset[g4decayChannelInstanceID]).parent_mass)
#define daughters_mass ((G4VDecayChannel::g4decayChannelSubInstanceManager.offset[g4decayChannelInstanceID]).daughters_mass)

#endif


class G4VDecayChannel
{
 // Class Description
 // This class is a abstract class to describe decay kinematics
 //

  public:

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.              
    //This new field is used as instance ID.                                        
    int g4decayChannelInstanceID;

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.              
    //This new field helps to use the class G4DecayChannelSubInstanceManager            
    //introduced above.                                                             
    static G4DecayChannelSubInstanceManager g4decayChannelSubInstanceManager;

    //Constructors 
      G4VDecayChannel(const G4String &aName, G4int Verbose = 1);
      G4VDecayChannel(const G4String  &aName, 
		     const G4String& theParentName,
		     G4double        theBR,
		     G4int           theNumberOfDaughters,
		     const G4String& theDaughterName1,
		     const G4String& theDaughterName2 = "",
		     const G4String& theDaughterName3 = "",
		     const G4String& theDaughterName4 = "" );

    //  Destructor
      virtual ~G4VDecayChannel();

  protected:
    //  default constructor
     G4VDecayChannel();
    //  copy constructor and assignment operatotr
     G4VDecayChannel(const G4VDecayChannel &);
     G4VDecayChannel & operator=(const G4VDecayChannel &);

  public:
    // equality operators
    G4int operator==(const G4VDecayChannel &right) const {return (this == &right);};
    G4int operator!=(const G4VDecayChannel &right) const {return (this != &right);};

    // less-than operator is defined for G4DecayTable
    G4int operator<(const G4VDecayChannel &right) const;

  public: // With Description
   virtual G4DecayProducts* DecayIt(G4double parentMass = -1.0) = 0;

  public: // With Description
     //get kinematics name
     const G4String&  GetKinematicsName() const;
     //get branching ratio
     G4double   GetBR() const;
     //get number of daughter particles
     G4int      GetNumberOfDaughters() const;     

     //get the pointer to the parent particle
     G4ParticleDefinition * GetParent();
     //get the pointer to a daughter particle 
     G4ParticleDefinition * GetDaughter(G4int anIndex);

     //get the angular momentum of the decay
     G4int GetAngularMomentum();
     //get the name of the parent particle
     const G4String& GetParentName() const;
     //get the name of a daughter particle
     const G4String& GetDaughterName(G4int anIndex) const;

     // get mass of parent
     G4double GetParentMass() const; 
     G4double GetDaughterMass(G4int anIndex) const; 

     //set the parent particle (by name or by pointer) 
     void SetParent(const G4ParticleDefinition * particle_type);
     void SetParent(const G4String &particle_name);
     //set branching ratio
     void SetBR(G4double value); 
     //set number of daughter particles
     void SetNumberOfDaughters(G4int value);     
     //set a daughter particle (by name or by pointer) 
     void SetDaughter(G4int anIndex, 
		      const G4ParticleDefinition * particle_type);
     void SetDaughter(G4int anIndex, 
		      const G4String &particle_name);

  protected: 
    // kinematics name
    G4String   kinematics_name;
    // branching ratio  [0.0 - 1.0]
    G4double   rbranch;
    // number of daughters
    G4int      numberOfDaughters;
    // parent particle
    G4String*  parent_name;
    //daughter particles
    G4String** daughters_name;

  protected: // With Description
    // celar daughters array
    void ClearDaughtersName();

  protected:
    // pointer to particle table
    G4ParticleTable*       particletable;

    // temporary buffers of pointers to G4ParticleDefinition
    // Change to thread private fields
    //    G4ParticleDefinition*  parentG4MTThreadPrivate;
    // Change to thread private fields
    //    G4ParticleDefinition** daughtersG4MTThreadPrivate;

    // parent mass
    // Change to thread private fields
    //    G4double               parent_mass;
    // Change to thread private fields
    //    G4double*              daughters_mass;

    // fill daughters array
    void FillDaughters();
    // fill parent
    void FillParent();

  public:  // With Description
    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel()  const;
    void  DumpInfo();

  private:
    const G4String& GetNoName() const;

  protected:  
    // controle flag for output message
    G4int verboseLevel;
    //  0: Silent
    //  1: Warning message
    //  2: More

    static const G4String   noName;
};

inline
 G4int G4VDecayChannel::operator<(const G4VDecayChannel &right) const
{
  return (this->rbranch < right.rbranch);
}

inline 
  G4ParticleDefinition* G4VDecayChannel::GetDaughter(G4int anIndex)
 { 
  //pointers to daughter particles are filled, if they are not set yet 
  if (daughtersG4MTThreadPrivate == 0) FillDaughters();

  //get the pointer to a daughter particle
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return daughtersG4MTThreadPrivate[anIndex];
  } else {
    if (verboseLevel>0)
      G4cout << "G4VDecayChannel::GetDaughter  index out of range "<<anIndex<<G4endl;
    return 0;
  }
}

inline
 const G4String& G4VDecayChannel::GetDaughterName(G4int anIndex) const
{
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return *daughters_name[anIndex];
  } else {
    if (verboseLevel>0){
      G4cout << "G4VDecayChannel::GetDaughterName ";
      G4cout << "index out of range " << anIndex << G4endl;
    }
    return GetNoName();
  }
}

inline
 G4double G4VDecayChannel::GetDaughterMass(G4int anIndex) const
{
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return daughters_mass[anIndex];
  } else {
    if (verboseLevel>0){
      G4cout << "G4VDecayChannel::GetDaughterMass ";
      G4cout << "index out of range " << anIndex << G4endl;
    }
    return 0.0;
  }
}

inline 
  G4ParticleDefinition* G4VDecayChannel::GetParent()
{ 
  //the pointer to the parent particle is filled, if it is not set yet 
   if (parentG4MTThreadPrivate == 0) FillParent();
  //get the pointer to the parent particle
  return parentG4MTThreadPrivate;
}

inline
 const G4String& G4VDecayChannel::GetParentName() const
{
  return *parent_name;
}

inline
 G4double G4VDecayChannel::GetParentMass() const
{
  return parent_mass;
}


inline
  void G4VDecayChannel::SetParent(const G4String &particle_name)
{
  if (parent_name != 0) delete parent_name;
  parent_name = new G4String(particle_name);
  parentG4MTThreadPrivate = 0;
}

inline
 G4int G4VDecayChannel::GetNumberOfDaughters() const 
{ 
  return  numberOfDaughters;
}

inline
 const G4String& G4VDecayChannel::GetKinematicsName() const { return kinematics_name; }

inline
 void  G4VDecayChannel::SetBR(G4double value){ rbranch = value; }

inline
 G4double G4VDecayChannel::GetBR() const { return rbranch; }

inline
 void  G4VDecayChannel::SetVerboseLevel(G4int value){ verboseLevel = value; }

inline
 G4int G4VDecayChannel::GetVerboseLevel() const { return verboseLevel; }



#endif






