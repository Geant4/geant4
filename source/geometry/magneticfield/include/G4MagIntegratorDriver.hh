// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagIntegratorDriver.hh,v 1.7 2000-11-01 15:15:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MagInt_Driver
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that 
// the error is within acceptable bounds.

// History:
// - Created. J.Apostolakis.

#ifndef G4MagInt_Driver_Def
#define G4MagInt_Driver_Def

#include "globals.hh"
#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"

class G4MagInt_Driver
{
   public:  // with description

     G4bool  AccurateAdvance(G4FieldTrack&  y_current,
		             G4double hstep,
			     G4double eps); // Requested y_err/hstep
       // Above drivers for integrator (Runge-Kutta) with stepsize control. 
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps. 
       // On output ystart is replaced by value at end of interval. 
       // The concept is similar to the odeint routine from NRC p.721-722.

     G4bool  QuickAdvance(G4FieldTrack& y_val,            // INOUT
			  const G4double     dydx[],  
			        G4double     hstep,       // IN 
				G4double&    dchord_step,
				G4double&    dyerr )  ;
        // QuickAdvance just tries one Step - it does not ensure accuracy.

     G4MagInt_Driver( G4double                hminimum, 
		      G4MagIntegratorStepper *pItsStepper,
                      G4int                   numberOfComponents=6);
     ~G4MagInt_Driver();
        // Constructor, destructor.

     inline G4double GetHmin() const;
     inline G4double Hmin() const;     // Obsolete
     inline G4double GetSafety() const;
     inline G4double GetPshrnk() const;
     inline G4double GetPgrow() const;
     inline G4double GetErrcon() const;
     inline void GetDerivatives( const G4FieldTrack y_curr,     // const, INput
			               G4double    dydx[]   );  //       OUTput
        // Accessors.

     inline void RenewStepperAndAdjust(G4MagIntegratorStepper *pItsStepper);
        // Sets a new stepper pItsStepper for this driver. Then it calls
	// ReSetParameters to reset its parameters accordingly.

     inline void ReSetParameters(G4double new_safety= 0.9 );
        //  i) sets the exponents (pgrow & pshrnk), 
        //     using the current Stepper's order, 
        // ii) sets the safety
        // ii) calculates "errcon" according to the above values.

     inline void SetSafety(G4double valS);
     inline void SetPshrnk(G4double valPs);
     inline void SetPgrow (G4double valPg);
     inline void SetErrcon(G4double valEc);
        // When setting safety or pgrow, errcon will be set to a 
        // compatible value.

     inline G4double ComputeAndSetErrcon();

     inline void SetChargeMomentumMass( G4double particleCharge,
		                        G4double MomentumXc,
		                        G4double Mass );
        // Change them in Equation. particleCharge is in e+ units.

     inline const G4MagIntegratorStepper* GetStepper() const;

     void  OneGoodStep(       G4double  ystart[], // Like old RKF45step()
			const G4double  dydx[],
			      G4double& x,
			      G4double htry,
			      G4double  eps,      //  memb variables ?
			      G4double& hdid,
			      G4double& hnext ) ;
        // This takes one Step that is as large as possible while 
        // satisfying the accuracy criterion of:
	// yerr < eps * |y_end-y_start|

     G4double ComputeNewStepSize( G4double  errMaxNorm,    // normalised error
			          G4double  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Do not limit the next step's size within a factor of the
	// current one.

     G4double ComputeNewStepSize_WithinLimits(
                          G4double  errMaxNorm,    // normalised error
			  G4double  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Limit the next step's size within a range around the current one.

     inline G4int    GetMaxNoSteps() const;
     inline void     SetMaxNoSteps( G4int val); 

   public:  // without description

     inline void SetHmin(G4double newval);

   protected:

     void WarnSmallStepSize( G4double hnext, G4double hstep, 
			     G4double h,     G4double xDone,
			     G4int noSteps);
     void WarnTooManyStep( G4double x1start, G4double x2end, G4double xCurrent);
     void WarnEndPointTooFar (G4double  endPointDist, 
			      G4double  hStepSize , 
			      G4double  epsilonRelative,
			      G4int     debugFlag);
        //  Issue warnings for undesirable situations

   private:

     G4MagInt_Driver(const G4MagInt_Driver&);
     G4MagInt_Driver& operator=(const G4MagInt_Driver&);
        // Private copy constructor and assignment operator.

   private:

     G4double hminimum_val;
        // Minimum Step allowed in a Step.

     G4int nvar;

     G4MagIntegratorStepper *pIntStepper;
     G4int   fMaxNoSteps;
     static const G4int  fMaxStepBase;  

     G4double safety;
     G4double pshrnk;   //  exponent for shrinking
     G4double pgrow;    //  exponent for growth
     G4double errcon;
        // Parameters used to grow and shrink trial stepsize.

     static const G4double max_stepping_increase;
     static const G4double max_stepping_decrease;
        // Maximum stepsize increase/decrease factors.
};

#include "G4MagIntegratorDriver.icc"

#endif /* G4MagInt_Driver_Def */
