// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronicInteraction.hh,v 1.1 1999-01-07 16:11:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Interaction  abstract base class
 // This class is the base class for the model classes.
 // It sorts out the energy-range for the models and provides
 // class utilities.
 // original by H.P. Wellisch
 // Modified by J.L.Chuma, TRIUMF, 21-Mar-1997
 // Last modified: 3-Apr-1997
 // Added units to energy initialization: J.L. Chuma  04-Apr-97
 // Modified by J.L.Chuma, 05-May-97  to Initialize theBlockedCounter
 // Modified by J.L.Chuma, 08-Jul-97 to implement the Nucleus changes
 
#ifndef G4HadronicInteraction_h
#define G4HadronicInteraction_h 1
 
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ReactionDynamics.hh"
#include "G4Material.hh"
#include "G4Nucleus.hh"
#include "G4Track.hh"
 
 class G4HadronicInteraction
 {
 public:
    
    G4HadronicInteraction() :
      verboseLevel(0), theMinEnergy(0.0*GeV), theMaxEnergy(25.0*GeV),
      theBlockedCounter(0), theMinCounter(0), theMaxCounter(0),
      theMinCounterElements(0), theMaxCounterElements(0),
      theBlockedCounterElements(0)
    { }
    
    virtual ~G4HadronicInteraction()
    { }
    
 private:
    
    inline G4HadronicInteraction(
     const G4HadronicInteraction &right )
    { *this = right; }
    
    inline const G4HadronicInteraction & operator=(
     const G4HadronicInteraction &right )
    { 
     if(this!=&right) G4Exception("unintended use of G4HadronicInteraction::operator=");
     return right;
    }
    
 public:
    
    inline G4bool operator==(
     const G4HadronicInteraction &right ) const
    { return ( this == (G4HadronicInteraction *) &right ); }
    
    inline G4bool operator!=(
     const G4HadronicInteraction &right ) const
    { return ( this != (G4HadronicInteraction *) &right ); }
    
    inline G4double GetMinEnergy() const
    { return theMinEnergy; }
    
    G4double GetMinEnergy( const G4Material *aMaterial,
                           const G4Element *anElement ) const;
    
    inline void SetMinEnergy( const G4double anEnergy )
    { theMinEnergy = anEnergy; }
    
    void SetMinEnergy( G4double anEnergy,
                       G4Element *anElement );
    
    void SetMinEnergy( G4double anEnergy,
                       G4Material *aMaterial );
    
    inline G4double GetMaxEnergy() const
    { return theMaxEnergy; }
    
    G4double GetMaxEnergy( const G4Material *aMaterial,
                           const G4Element *anElement ) const;
    
    inline void SetMaxEnergy( const G4double anEnergy )
    { theMaxEnergy = anEnergy; }
    
    void SetMaxEnergy( G4double anEnergy,
                       G4Element *anElement );
    
    void SetMaxEnergy( G4double anEnergy,
                       G4Material *aMaterial );
  
    inline const G4HadronicInteraction *GetMyPointer() const
    { return this; }

    inline G4int GetVerboseLevel() const
    { return verboseLevel; }

    inline void SetVerboseLevel( G4int value )
    { verboseLevel = value; }

    virtual G4VParticleChange *ApplyYourself(
     const G4Track &aTrack, G4Nucleus & targetNucleus ) = 0;

    void DeActivateFor( G4Material *aMaterial );

    void DeActivateFor( G4Element *anElement ); 

    G4bool IsBlocked( const G4Material *aMaterial ) const;

    G4bool IsBlocked( const G4Element *anElement) const;
    
 protected:
    
    G4ParticleChange theParticleChange;
    // the G4VParticleChange object which is modified and returned
    // by address by the ApplyYourself method,
    // (instead of aParticleChange as found in G4VProcess)
    
    G4int verboseLevel;
    // control flag for output messages
    // 0: silent
    // 1: warning messages
    // 2: more
    // (instead of verboseLevel as found in G4VProcess)
    
    G4ReactionDynamics theReactionDynamics;
    
    // these two have global validity
    // units are assumed to be MeV
    
    G4double theMinEnergy;
    G4double theMaxEnergy;
    
 private:
    
    enum { MAX_LIST_SIZE = 500 };
    
    // the following allow for restrictions/additions for specific materials
    
    G4double theMinEnergyList[ MAX_LIST_SIZE ];
    G4Material *theMinMaterials[ MAX_LIST_SIZE ];
    G4int theMinCounter;
    
    G4double theMaxEnergyList[ MAX_LIST_SIZE ];
    G4Material *theMaxMaterials[ MAX_LIST_SIZE ];
    G4int theMaxCounter;

    G4Material *theBlockedList[ MAX_LIST_SIZE ];
    G4int theBlockedCounter;

    // the following allow for restrictions/additions for specific elements

    G4double theMinEnergyListElements[ MAX_LIST_SIZE ];
    G4Element *theMinElements[ MAX_LIST_SIZE ];
    G4int theMinCounterElements;
    
    G4double theMaxEnergyListElements[ MAX_LIST_SIZE ];
    G4Element *theMaxElements[ MAX_LIST_SIZE ];
    G4int theMaxCounterElements;

    G4Element *theBlockedListElements[ MAX_LIST_SIZE ];
    G4int theBlockedCounterElements;
 };
 
#endif
