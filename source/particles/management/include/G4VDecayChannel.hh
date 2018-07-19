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
// $Id: G4VDecayChannel.hh 105720 2017-08-16 12:38:10Z gcosmo $
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

#include <cmath>

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

class    G4ParticleDefinition;
class    G4DecayProducts;
class    G4ParticleTable;

class G4VDecayChannel
{
  // Class Description:
  //
  // Abstract class to describe decay kinematics

  public:  // with description

    // Constructors 
    G4VDecayChannel(const G4String &aName, G4int Verbose = 1);
    G4VDecayChannel(const G4String  &aName, 
		    const G4String& theParentName,
		    G4double        theBR,
		    G4int           theNumberOfDaughters,
		    const G4String& theDaughterName1,
		    const G4String& theDaughterName2 = "",
		    const G4String& theDaughterName3 = "",
		    const G4String& theDaughterName4 = "" );

    // Destructor
    virtual ~G4VDecayChannel();

    // equality operators
    G4int operator==(const G4VDecayChannel &right) const {return (this == &right);}
    G4int operator!=(const G4VDecayChannel &right) const {return (this != &right);}

    // less-than operator is defined for G4DecayTable
    G4int operator<(const G4VDecayChannel &right) const;

    virtual G4DecayProducts* DecayIt(G4double parentMass = -1.0) = 0;

    // get kinematics name
    const G4String&  GetKinematicsName() const;
    // get branching ratio
    G4double   GetBR() const;
    // get number of daughter particles
    G4int      GetNumberOfDaughters() const;     

    // get the pointer to the parent particle
    G4ParticleDefinition * GetParent();
    // get the pointer to a daughter particle 
    G4ParticleDefinition * GetDaughter(G4int anIndex);

    // get the angular momentum of the decay
    G4int GetAngularMomentum();
    // get the name of the parent particle
    const G4String& GetParentName() const;
    //get the name of a daughter particle
    const G4String& GetDaughterName(G4int anIndex) const;

    // get mass of parent
    G4double GetParentMass() const; 
    G4double GetDaughterMass(G4int anIndex) const; 

    // set the parent particle (by name or by pointer) 
    void SetParent(const G4ParticleDefinition * particle_type);
    void SetParent(const G4String &particle_name);
    // set branching ratio
    void SetBR(G4double value); 
    // set number of daughter particles
    void SetNumberOfDaughters(G4int value);     
    // set a daughter particle (by name or by pointer) 
    void SetDaughter(G4int anIndex, 
                     const G4ParticleDefinition * particle_type);
    void SetDaughter(G4int anIndex, 
                     const G4String &particle_name);

    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel()  const;
    void  DumpInfo();
 
    G4double  GetRangeMass() const;
    void      SetRangeMass(G4double val);
    virtual G4bool    IsOKWithParentMass(G4double parentMass); 

    void  SetPolarization(const G4ThreeVector&);
    const G4ThreeVector& GetPolarization() const;

  protected: // with description

    // clear daughters array
    void ClearDaughtersName();

  protected:
    inline void CheckAndFillDaughters();
    inline void CheckAndFillParent();

  private:
    // fill daughters array
    void FillDaughters();
    // fill parent
    void FillParent();

  protected: // without description

    // default constructor
    G4VDecayChannel();

    // copy constructor and assignment operator
    G4VDecayChannel(const G4VDecayChannel &);
    G4VDecayChannel & operator=(const G4VDecayChannel &);

  private:

    const G4String& GetNoName() const;

  protected:
     G4double DynamicalMass(G4double massPDG, G4double width, G4double maxDev = +1.) const;

  protected:

    // kinematics name
    G4String   kinematics_name;
    // branching ratio  [0.0 - 1.0]
    G4double   rbranch;
    // number of daughters
    G4int      numberOfDaughters;
    // parent particle
    G4String*  parent_name;
    // daughter particles
    G4String** daughters_name;
 
    //  range of mass allowed in decay
    G4double   rangeMass;

    //  polarization of the parent particle
    G4ThreeVector parent_polarization;
   
    // pointer to particle table
    G4ParticleTable*       particletable;

    // control flag for output message
    G4int verboseLevel;
    //  0: Silent
    //  1: Warning message
    //  2: More

    static const G4String noName;

    G4ParticleDefinition*  G4MT_parent;
    G4ParticleDefinition** G4MT_daughters;
    G4double               G4MT_parent_mass;
    G4double*              G4MT_daughters_mass;
    G4double*              G4MT_daughters_width;
    G4Mutex daughtersMutex;
    G4Mutex parentMutex;
};

// ------------------------------------------------------------
// Inline methods
// ------------------------------------------------------------

inline
 G4int G4VDecayChannel::operator<(const G4VDecayChannel &right) const
{
  return (this->rbranch < right.rbranch);
}

inline 
  G4ParticleDefinition* G4VDecayChannel::GetDaughter(G4int anIndex)
 { 
  // pointers to daughter particles are filled, if they are not set yet 
  CheckAndFillDaughters();

  // get the pointer to a daughter particle
  if ( (anIndex>=0) && (anIndex<numberOfDaughters) ) {
    return G4MT_daughters[anIndex];
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
    return G4MT_daughters_mass[anIndex];
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
  // the pointer to the parent particle is filled, if it is not set yet 
  CheckAndFillParent();
  // get the pointer to the parent particle
  return G4MT_parent;
}

inline
 const G4String& G4VDecayChannel::GetParentName() const
{
  return *parent_name;
}

inline
 G4double G4VDecayChannel::GetParentMass() const
{
  return G4MT_parent_mass;
}

inline
  void G4VDecayChannel::SetParent(const G4String &particle_name)
{
  if (parent_name != 0) delete parent_name;
  parent_name = new G4String(particle_name);
  G4MT_parent = 0;
}

inline
 G4int G4VDecayChannel::GetNumberOfDaughters() const 
{ 
  return  numberOfDaughters;
}

inline
 const G4String& G4VDecayChannel::GetKinematicsName() const { return kinematics_name; }

inline
 G4double G4VDecayChannel::GetBR() const { return rbranch; }

inline
 void  G4VDecayChannel::SetVerboseLevel(G4int value){ verboseLevel = value; }

inline
 G4int G4VDecayChannel::GetVerboseLevel() const { return verboseLevel; }

inline   
 G4double  G4VDecayChannel::GetRangeMass() const { return rangeMass; }

inline   
  void  G4VDecayChannel::SetRangeMass(G4double val){ if(val>=0.) rangeMass=val; }

inline 
  void G4VDecayChannel::SetPolarization(const G4ThreeVector& polar)
{
  parent_polarization = polar;
}

inline 
  const G4ThreeVector& G4VDecayChannel::GetPolarization() const
{
  return parent_polarization;
}

inline
void G4VDecayChannel::CheckAndFillDaughters() {
  G4AutoLock l(&daughtersMutex);
  if ( G4MT_daughters ==0 ) {
      l.unlock();
      FillDaughters();
  }
}

inline void G4VDecayChannel::CheckAndFillParent() {
  G4AutoLock l(&parentMutex);
  if ( G4MT_parent == 0 ) {
      l.unlock();
      FillParent();
  }
}
#endif
