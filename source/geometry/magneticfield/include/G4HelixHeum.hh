// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixHeum.hh,v 1.4 2000-04-27 09:14:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4HelixHeum
//
// Class description:
//
// Simple Heum:
//     x_1 = x_0 +
//           h * 1/4 * dx(t0,x0)  +
//               3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 

// History:
// - Created. W.Wander <wwc@mit.edu>, 03/11/98

#ifndef G4HELIXHEUM_HH
#define G4HELIXHEUM_HH

#include "G4MagHelicalStepper.hh"

class G4HelixHeum : public G4MagHelicalStepper
{

  public:  // with description

    G4HelixHeum(G4Mag_EqRhs *EqRhs)
      : G4MagHelicalStepper(EqRhs) {;}

    ~G4HelixHeum() {;}
  
    void DumbStepper( const G4double y[],
		      G4ThreeVector  Bfld,
		      G4double       h,
		      G4double       yout[]);

  public: // without description
  
    G4int IntegratorOrder() { return 2; }
};

#endif /* G4HELIXHEUM_HH */
