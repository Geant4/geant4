// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixImplicitEuler.hh,v 1.4 2000-04-27 09:14:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4HelixImplicitEuler
//
// Class description:
//
//  Helix Implicit Euler:
//        x_1 = x_0 + 1/2 * ( helix(h,t_0,x_0)
//                          + helix(h,t_0+h,x_0+helix(h,t0,x0) ) )

// History:
// - Created. W.Wander <wwc@mit.edu>, 03/11/98

#ifndef G4HELIXIMPLICITEULER_HH
#define G4HELIXIMPLICITEULER_HH

#include "G4MagHelicalStepper.hh"

class G4HelixImplicitEuler : public G4MagHelicalStepper
{

  public:  // with description

    G4HelixImplicitEuler(G4Mag_EqRhs *EqRhs)
      : G4MagHelicalStepper(EqRhs) {;}

    ~G4HelixImplicitEuler() {;}
  
    void DumbStepper( const G4double y[],
		      G4ThreeVector  Bfld,
		      G4double       h,
		      G4double       yout[]);
  
  public:  // without description

    G4int IntegratorOrder() { return 2; }
};

#endif /* G4HELIXIMPLICITEULER_HH */
