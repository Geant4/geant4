// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixExplicitEuler.hh,v 1.3 2000-04-12 18:28:50 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// W. Wander <wwc@mit.edu> 12/09/97

#ifndef G4HELIXEXPLICITEULER_HH
#define G4HELIXEXPLICITEULER_HH
#include "G4MagHelicalStepper.hh"

class G4HelixExplicitEuler: public G4MagHelicalStepper
{

public:
  G4HelixExplicitEuler(G4Mag_EqRhs *EqRhs): G4MagHelicalStepper(EqRhs){};
  ~G4HelixExplicitEuler(){};
  
  void  DumbStepper(  const G4double y[],
		      G4ThreeVector  Bfld,
		      G4double       h,
		      G4double       yout[]);

  // DELETED  RightHandSide( ) !!!!  
  // Replace by MagFieldEvaluate( const G4double y[], G4double B[] )   
  //  in G4HelicalStepper
  
  G4int     IntegratorOrder() { return 1; };
};

#endif /* G4EXPLICITEULER_HH */
