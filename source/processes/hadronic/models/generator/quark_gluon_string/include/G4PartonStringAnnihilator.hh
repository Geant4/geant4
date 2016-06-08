// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PartonStringAnnihilator.hh,v 1.3 1999/12/15 14:52:42 gunter Exp $
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 9-Oct-1998
// -----------------------------------------------------------------------------
#ifndef G4PartonStringAnnihilator_h
#define G4PartonStringAnnihilator_h 1

#include "G4KineticTrack.hh"
#include "G4ExcitedString.hh"
#include "G4MesonSplitter.hh"
#include "G4BaryonSplitter.hh"


class G4PartonStringAnnihilator
    {   
public:
    G4PartonStringAnnihilator();
    
public:
   
    G4double GetCrossSection(G4KineticTrack& aTarget, G4KineticTrack& aProjectile);
    G4ExcitedString* GetString(G4KineticTrack& Target, G4KineticTrack& Projectile);

    G4bool SplitUpBarion(G4int Encoding, G4int* q_or_qqbar, G4int* qbar_or_qq); //!
 
    G4bool isMeson(G4int Encoding); 
    G4bool isBarion(G4int Encoding);
    G4bool isAntiParticle(G4int Encoding);

private:    
    G4bool GetStringEnds(G4int Projectile, G4int Target, G4int* Left, G4int* Right);
    G4bool FindDiquark(G4int Encoding, G4int Quark, G4int* Diquark);
    G4bool FindQuark(G4int Encoding, G4int Diquark, G4int* Quark);

private:

    G4double widthOfPtSquare;
    G4MesonSplitter theMesonSplitter;
    G4BaryonSplitter theBaryonSplitter;
    };

// *************************************************************************************************
#endif
