// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PartonStringAnnihilator.hh,v 1.1.4.1 1999/12/07 20:51:51 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $Maxim Komogorov
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 9-Oct-1998
// -----------------------------------------------------------------------------
#ifndef G4PartonStringAnnihilator_h
#define G4PartonStringAnnihilator_h 1

#include "G4KineticTrackVector.hh"
#include "G4ExcitedString.hh"
#include "G4MesonSplitter.hh"

// **************************************************************************************************

// Barion consists from Quark and Diquark with given spin-isospin state.
// Different barion states are defined by Probability

struct G4ACSParameters
    {
    G4int ProjectileEncoding;
    G4int TargetEncoding;
    G4double X;
    G4double Y;
    G4double Eta;
    G4double Epsilon;
    };

// **************************************************************************************************

class G4PartonStringAnnihilator
    {   
public:
    G4PartonStringAnnihilator();
   ~G4PartonStringAnnihilator();
    static G4ACSParameters ACSParametersTable[];    
    
public:
   
    G4bool IsAnnihilation(G4KineticTrackVector& aTarget, G4KineticTrackVector& aProjectile, G4double* CrossSection);
    G4ExcitedString* Annihilator(G4KineticTrack& Target, G4KineticTrack& Projectile);

private:    
//    void   GetValenceQuarkFlavors(G4int PDGcode, G4int& aEnd, G4int& bEnd);
    G4bool AnnihilatorEncoding(G4int Projectile, G4int Target, G4int* Left, G4int* Right);
    G4bool FindDiquark(G4int Encoding, G4int Quark, G4int* Diquark);
    G4bool FindQuark(G4int Encoding, G4int Diquark, G4int* Quark);
public:
    G4bool SplitUpBarion(G4int Encoding, G4int* q_or_qqbar, G4int* qbar_or_qq); //!
 
    G4bool isMeson(G4int Encoding); 
    G4bool isBarion(G4int Encoding);
    G4bool isAntiParticle(G4int Encoding);
private:

    G4double widthOfPtSquare;
    G4MesonSplitter theMesonSplitter;
    };

// *************************************************************************************************
#endif
