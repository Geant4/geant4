// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagIntegratorDriver.hh,v 1.4 1999-12-15 14:49:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  Provides Driver that talks to Integrator Stepper, and insures that 
//  the error is within acceptable bounds.

#ifndef G4MagInt_Driver_Def
#define G4MagInt_Driver_Def

#include "globals.hh"
#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"

class G4MagInt_Driver
{
   public:
      G4bool  AccurateAdvance(    G4FieldTrack&  y_current,
			    const G4double      hstep,
			    const G4double      eps); // Requested y_err/hstep

      //
      // Above drivers for integrator (Runge-Kutta) with stepsize control. 
      //  Integrates ODE starting values y_current
      // from current s (s=s0) to s=s0+h with accuracy eps. 
      // On output ystart is replaced by value at end of interval. 
      // The concept is similar to the odeint routine from NRC p.721-722 .

      // QuickAdvance just tries one Step - it does not ensure accuracy
      G4bool  QuickAdvance(       G4FieldTrack& y_val,   // INOUT
			    const G4double     dydx[],  
			          G4double     hstep,       // IN 
				  G4double&    dchord_step,
				  G4double&    dyerr )  ;

      // Constructor, destructor
      // 
      G4MagInt_Driver( G4double                hminimum, 
		       G4MagIntegratorStepper *pItsStepper,
                       G4int                   numberOfComponents=6);
     ~G4MagInt_Driver(){}

      // Access functions
      // ----------------
      //
      G4double GetHmin(){ return hminimum_val;} 
      G4double Hmin()   { return hminimum_val;}     // Obsolete
      G4double GetSafety(){ return safety; }
      G4double GetPshrnk(){ return pshrnk;} 
      G4double GetPgrow(){ return pgrow;} 
      G4double GetErrcon(){ return errcon;}

      void   GetDerivatives( const G4FieldTrack y_curr,     // const, INput
			           G4double    dydx[]   );  //       OUTput

      // Set    functions
      // ----------------
      //
      G4double SetHmin(G4double newval){ return hminimum_val;} 
      //
      //  The following function sets a new stepper pItsStepper for
      //   this driver, and then calls ReSetParameters to reset its 
      //   parameters accordingly.
      //
      void RenewStepperAndAdjust(G4MagIntegratorStepper *pItsStepper);
      //
      //  
      //   ReSetParameters does the following: 
      //    i) sets the exponents (pgrow & pshrnk), 
      //                    using the current Stepper's order, 
      //   ii) sets the safety
      //   ii) calculates "errcon" according to the above values.
      //
      void ReSetParameters(G4double new_safety= 0.9 );
      //
      // When setting safety or pgrow, errcon will be set to a 
      //   compatible value 
      //
      void SetSafety(G4double valS);
      void SetPshrnk(G4double valPs);
      void SetPgrow( G4double valPg);
      void SetErrcon(G4double valEc);
      //
      G4double ComputeAndSetErrcon();


      void   SetChargeMomentumMass(                // Change them in Equation
		     const G4double particleCharge,    // in e+ units
		     const G4double MomentumXc,
		     const G4double Mass );

      G4MagIntegratorStepper* GetStepper();

      //  This takes one Step that is as large as possible while 
      //  satisfying the accuracy criterion of 
      //                 yerr < eps * |y_end-y_start|
      //
      void  OneGoodStep(       G4double  ystart[], // Like old RKF45step()
			 const G4double  dydx[],
			       G4double& x,
			 const G4double htry,
			 const G4double  eps,      //  memb variables ?
			       G4double& hdid,
			       G4double& hnext ) ;

      //  Taking the last step's normalised error, calculate
      //    a step size for the next step.
      //  Do not limit the next step's size within a factor of the current one.
      G4double ComputeNewStepSize( 
			  G4double  errMaxNorm,    // normalised error
			  G4double  hstepCurrent); // current step size
                                          
      //  Taking the last step's normalised error, calculate
      //    a step size for the next step.
      //  Limit the next step's size within a range around the current one.
      G4double ComputeNewStepSize_WithinLimits( 
			  G4double  errMaxNorm,    // normalised error
			  G4double  hstepCurrent); // current step size

      G4int    GetMaxNoSteps();
      void     SetMaxNoSteps( G4int val); 

protected:
      //  Issue warnings for undesirable situations
      void WarnSmallStepSize( G4double hnext, G4double hstep, 
			      G4double h,     G4double xDone,
			      G4int noSteps);
      void WarnTooManyStep( G4double x1start, G4double x2end, G4double xCurrent);
      void WarnEndPointTooFar (G4double  endPointDist, 
			       G4double  hStepSize , 
			       G4double  epsilonRelative,
			       G4int     debugFlag);

private:
      G4double hminimum_val;  //  Minimum Step allowed in a Step

      const G4int nvar;

      G4MagIntegratorStepper *pIntStepper;
      G4int   fMaxNoSteps;
      static const G4int  fMaxStepBase;  

      //  Parameters used to grow and shrink trial stepsize
      G4double safety;
      G4double pshrnk;   //  exponent for shrinking
      G4double pgrow;    //  exponent for growth
      G4double errcon;
      //  maximum stepsize increase/decrease factors
      const static G4double max_stepping_increase;
      const static G4double max_stepping_decrease;
};

#include "G4MagIntegratorDriver.icc"

#endif /* G4MagInt_Driver_Def */
