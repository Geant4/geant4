// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyRangeManager.hh,v 1.2 1999-12-15 14:52:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Energy Range Manager
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Last modified: 24-Mar-1997
 
#ifndef G4EnergyRangeManager_h
#define G4EnergyRangeManager_h 1
 
#include "G4HadronicInteraction.hh"
 
 class G4EnergyRangeManager 
 {
     
 public:
    
    G4EnergyRangeManager()
    { theHadronicInteractionCounter = 0; }
    
    ~G4EnergyRangeManager()
    { }
    
    G4EnergyRangeManager( const G4EnergyRangeManager &right );
    
    G4EnergyRangeManager & operator=( const G4EnergyRangeManager &right );
    
 public:
    
    inline G4bool operator==( const G4EnergyRangeManager &right ) const
    { return ( this == (G4EnergyRangeManager *) &right ); }
    
    inline G4bool operator!=( const G4EnergyRangeManager &right ) const
    { return ( this != (G4EnergyRangeManager *) &right ); }
    
    void RegisterMe( G4HadronicInteraction *a );
    
    G4HadronicInteraction *GetHadronicInteraction(
     const G4double kineticEnergy,
     const G4Material *aMaterial,
     const G4Element *anElement ) const;
    
 private:
    
    inline G4int GetHadronicInteractionCounter() const
    { return theHadronicInteractionCounter; }
    
 private:
     
    enum { MAX_NUMBER_OF_MODELS = 100 };
    
    G4HadronicInteraction *
     theHadronicInteraction[ MAX_NUMBER_OF_MODELS ];
    
    G4int theHadronicInteractionCounter;
    
 };

#endif
 
