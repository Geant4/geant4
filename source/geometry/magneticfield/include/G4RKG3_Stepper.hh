// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RKG3_Stepper.hh,v 1.1 1999-01-07 16:07:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// J.Apostolakis, V.Grichine 30.01.97

#include "G4MagIntegratorStepper.hh"
#include "G4ThreeVector.hh"

class G4RKG3_Stepper : public G4MagIntegratorStepper
{
public:
        G4RKG3_Stepper(G4Mag_EqRhs *EqRhs): G4MagIntegratorStepper(EqRhs){};
       ~G4RKG3_Stepper(){};

	//  The method it must provide, even if less efficiently
	
	void
	Stepper( const G4double yIn[],
		 const G4double dydx[],
		 const G4double h,
		       G4double yOut[],
		       G4double yErr[]  );
		                                 //  G4double& beta2) const

        G4double  DistChord() const ;
 
         // Additional "optimised" methods: 

	// Integrator RK Stepper from G3 with only two field evaluation per 
	// Step. It is used in propagation initial Step by small substeps
	// after solution error and delta geometry considerations. 
	// B[3] is magnetic field which is passed from substep to substep.

        void StepNoErr(  const G4double tIn[7],
			 const G4double dydx[7],
			 const G4double Step,
			       G4double tOut[7],
			       G4double B[3] );

        void StepWithEst(const G4double  tIn[7],
			 const G4double dydx[7],
			 const G4double Step,
			       G4double tOut[7],
			                            //  G4double tError[6],
			       G4double& alpha2,    // to delete ?
			       G4double& beta2,
			 const G4double B1[3],
			       G4double B3[3] );

        G4int      IntegratorOrder() { return 4; };
  
protected:
  //	void Field( const  double Point[3],
  //				   double Bfield[3] )  const
  //			   { EqRhs-> GetFieldValue( Point, Bfield ) ; }
private:
          G4ThreeVector fyInitial,
	                fyMidPoint,
			fyFinal     ;
};

