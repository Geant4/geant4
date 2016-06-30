//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4FSALIntegrationDriver.hh 97387 2016-06-02 10:03:42Z gcosmo $
//
//
// class G4FSALIntegrationDriver
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that 
// the error is within acceptable bounds.

// History:
// - Created. J.Apostolakis.
// --------------------------------------------------------------------

#ifndef G4FSALIntegrationDriver_Def
#define G4FSALIntegrationDriver_Def

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4VFSALIntegrationStepper.hh"

class G4FSALIntegrationDriver
{
   public:  // with description

     G4bool  AccurateAdvance(G4FieldTrack&  y_current,
                             G4double hstep,
                             G4double eps,            // Requested y_err/hstep
                             G4double hinitial=0.0);  // Suggested 1st interval
       // Above drivers for integrator (Runge-Kutta) with stepsize control. 
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps. 
       // On output ystart is replaced by value at end of interval. 
       // The concept is similar to the odeint routine from NRC p.721-722.

     G4bool  QuickAdvance(      G4FieldTrack& y_val,      // INOUT
                           		G4double     dydx[],
                                G4double     hstep,       // IN 
                                G4double&    dchord_step,
                                G4double&    dyerr )  ;
        // QuickAdvance just tries one Step - it does not ensure accuracy.

     G4bool  QuickAdvance(      G4FieldTrack& y_posvel,        // INOUT
                                G4double      dydx[],
                                G4double      hstep,           // IN
                                G4double&     dchord_step,
                                G4double&     dyerr_pos_sq,
                                G4double&     dyerr_mom_rel_sq ) ;
       // New QuickAdvance that also just tries one Step
       //    (so also does not ensure accuracy)
       //    but does return the errors in  position and
       //        momentum (normalised: Delta_Integration(p^2)/(p^2) )

     G4FSALIntegrationDriver( G4double                hminimum, 
                        G4VFSALIntegrationStepper *pItsStepper,
                        G4int                   numberOfComponents=6,
                        G4int                   statisticsVerbosity=1);
     ~G4FSALIntegrationDriver();
        // Constructor, destructor.

     inline G4double GetHmin() const;
     inline G4double Hmin() const;     // Obsolete
     inline G4double GetSafety() const;
     inline G4double GetPshrnk() const;
     inline G4double GetPgrow() const;
     inline G4double GetErrcon() const;
     inline G4int GetNoTotalSteps() const;	//Only for debug purposes
     inline void GetDerivatives( const G4FieldTrack &y_curr,     // const, INput
                                       G4double    dydx[]   );  //       OUTput
        // Accessors.

     inline void RenewStepperAndAdjust(G4VFSALIntegrationStepper *pItsStepper);
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

     inline const G4VFSALIntegrationStepper* GetStepper() const;
     inline G4VFSALIntegrationStepper* GetStepper();

     void  OneGoodStep(       G4double  ystart[], // Like old RKF45step()
                        	  G4double  dydx[],
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
        //  Modify and Get the Maximum number of Steps that can be
        //   taken for the integration of a single segment -
        //   (ie a single call to AccurateAdvance).
    
//---------------------------------------------------------------------
//The following has been introduced by [hackabot] for testing purposes only
    inline G4int GetTotalNoStepperCalls() const;
//---------------------------------------------------------------------

   public:  // without description

     inline void SetHmin(G4double newval);
     inline void SetVerboseLevel(G4int newLevel); 
     inline G4double GetVerboseLevel() const;

     inline G4double GetSmallestFraction() const; 
     void     SetSmallestFraction( G4double val );
    
    
    

   protected:  // without description
     void WarnSmallStepSize( G4double hnext, G4double hstep, 
                             G4double h,     G4double xDone,
                             G4int noSteps);
     void WarnTooManyStep( G4double x1start, G4double x2end, G4double xCurrent);
     void WarnEndPointTooFar (G4double  endPointDist, 
                              G4double  hStepSize , 
                              G4double  epsilonRelative,
                              G4int     debugFlag);
        //  Issue warnings for undesirable situations

     void PrintStatus(  const G4double*      StartArr,
                              G4double       xstart,
                        const G4double*      CurrentArr, 
                              G4double       xcurrent, 
                              G4double       requestStep, 
                              G4int          subStepNo );
     void PrintStatus(  const G4FieldTrack&  StartFT,
                        const G4FieldTrack&  CurrentFT, 
                              G4double       requestStep, 
                              G4int          subStepNo );
     void PrintStat_Aux( const G4FieldTrack& aFieldTrack,
                               G4double      requestStep, 
                               G4double      actualStep,
                               G4int         subStepNo,
                               G4double      subStepSize,
                               G4double      dotVelocities );       
       //  Verbose output for debugging

     void PrintStatisticsReport() ;
       //  Report on the number of steps, maximum errors etc.

#ifdef QUICK_ADV_TWO
     G4bool QuickAdvance(      G4double     yarrin[],     // In
                         const G4double     dydx[],  
                               G4double     hstep,        
                               G4double     yarrout[],    // Out
                               G4double&    dchord_step,  // Out
                               G4double&    dyerr );      // in length
#endif

   private:

     G4FSALIntegrationDriver(const G4FSALIntegrationDriver&);
     G4FSALIntegrationDriver& operator=(const G4FSALIntegrationDriver&);
        // Private copy constructor and assignment operator.

   private:

     // ---------------------------------------------------------------
     //  INVARIANTS 

     G4double  fMinimumStep;
        // Minimum Step allowed in a Step (in absolute units)
     G4double  fSmallestFraction;      //   Expected range 1e-12 to 5e-15;  
        // Smallest fraction of (existing) curve length - in relative units
        //  below this fraction the current step will be the last 

     const G4int  fNoIntegrationVariables;  // Number of Variables in integration
     const G4int  fMinNoVars;               // Minimum number for FieldTrack
     const G4int  fNoVars;                  // Full number of variable

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

     G4int    fStatisticsVerboseLevel;

     // ---------------------------------------------------------------
     // DEPENDENT Objects
     G4VFSALIntegrationStepper *pIntStepper;

     // ---------------------------------------------------------------
     //  STATE

     G4int  fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps;
     G4double fDyerr_max, fDyerr_mx2;
     G4double fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
     G4double fSumH_sm, fSumH_lg; 
        // Step Statistics 

     G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
        // Could be varied during tracking - to help identify issues
    
    //For Test Purposes :-
    G4int TotalNoStepperCalls;

};

#include "G4FSALIntegrationDriver.icc"

#endif /* G4FSALIntegrationDriver_Def */
