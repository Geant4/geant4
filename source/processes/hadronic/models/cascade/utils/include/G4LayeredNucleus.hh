// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// original class G4Nucleus by H.P. Wellisch
//   modified by J.L. Chuma, TRIUMF, 19-Nov-1996
//   last modified: 27-Mar-1997
//   Chr. Volcker, 10-Nov-1997: new methods and class variables.
//   M.G. Pia, 2 Oct 1998: modified GetFermiMomentum (original design was
//                       the source of memory leaks)
// This class G4LayeredNucleus by T. Lampen 14.6.2000

 
#ifndef G4LayeredNucleus_h
#define G4LayeredNucleus_h 1
 
#include "G4Nucleus.hh"
 
class G4LayeredNucleus  : public G4Nucleus
{
public:
  
  G4LayeredNucleus() : G4Nucleus()
  { 
  }
  
  G4LayeredNucleus( const G4double A, const G4double Z ) : G4Nucleus ( A,  Z )
  {
  }

  G4LayeredNucleus( const G4Material *aMaterial ) : G4Nucleus ( aMaterial )
  {
  }

  G4ThreeVector GetMomentum();
  void SetMomentum( G4ThreeVector mom);

private:
    
  G4ThreeVector momentumVector;

};
 
#endif
 
