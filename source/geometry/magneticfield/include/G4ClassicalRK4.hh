// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ClassicalRK4.hh,v 1.2 1999-12-15 14:49:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// J.Apostolakis, V.Grichine 30.1.97
// changed: W.Wander <wwc@mit.edu> 12/09/97: Moved into MagErrorStepper

#include "G4MagErrorStepper.hh"
#include "G4ThreeVector.hh"

class G4ClassicalRK4 : public G4MagErrorStepper 
{
  public:

    G4ClassicalRK4(G4Mag_EqRhs *EqRhs, G4int numberOfVariables = 6) ;

   ~G4ClassicalRK4() ;
   
           void StepWithEst( const G4double  yIn[],
			     const G4double  dydx[],
			     const G4double  h,
			 	   G4double  yOut[],
                                   G4double& alpha2,
                                   G4double& beta2,
			     const G4double B1[],
			           G4double B2[]         ) ;

  
            G4int        IntegratorOrder() { return 4; };
  

   // A stepper that does not know about errors.
   //  It is used by the MagErrorStepper stepper.
   
            void DumbStepper( const G4double  yIn[],
			      const G4double  dydx[],
			      const G4double  h,
				    G4double  yOut[]) ;

   //  Could make above G4SixPoint to keep tangents too ...?

private:

  G4int fNumberOfVariables ; // is set default to 6 in constructor

  G4double *dydxm, *dydxt, *yt; // scratch space - not state 
};
