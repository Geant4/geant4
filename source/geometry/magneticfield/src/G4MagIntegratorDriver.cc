//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MagIntegratorDriver.cc,v 1.27 2002-06-11 08:15:52 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Implementation for class G4MagInt_Driver
//   Tracking in space dependent magnetic field
//
// History of major changes:
//  8 Nov 01  J. Apostolakis:   Respect minimum step in AccurateAdvance
// 27 Jul 99  J. Apostolakis:   Ensured that AccurateAdvance does not loop 
//                                due to very small eps & step size (precision)
// 28 Jan 98  W. Wander:        Added ability for low order integrators
//  7 Oct 96  V. Grichine       First version
#include <math.h>
#include "G4ios.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4FieldTrack.hh"
#include "geomdefs.hh"         //  for kCarTolerance
#include "g4std/iomanip"

//  Stepsize can increase by no more than 5.0
//           and decrease by no more than 1/10. = 0.1
//
const G4double G4MagInt_Driver::max_stepping_increase = 5.0;
const G4double G4MagInt_Driver::max_stepping_decrease = 0.1;

//  The (default) maximum number of steps is Base divided by the order of Stepper
//
const G4int  G4MagInt_Driver::fMaxStepBase = 250;  // Was 5000

//  Constructor
//
G4MagInt_Driver::G4MagInt_Driver( G4double                hminimum, 
				  G4MagIntegratorStepper *pItsStepper,
				  G4int                   numComponents)
  : nvar(numComponents), 
    fVerboseLevel(0)
{  
      RenewStepperAndAdjust( pItsStepper );
      hminimum_val= hminimum;
      fMaxNoSteps = fMaxStepBase / pIntStepper->IntegratorOrder();
#ifdef G4DEBUG_FIELD
      fVerboseLevel=2;
#endif
}

//  Destructor
//
G4MagInt_Driver::~G4MagInt_Driver()
{
}

// To add much printing for debugging purposes, uncomment this:
// #define  G4DEBUG_FIELD 1    

G4bool
G4MagInt_Driver::AccurateAdvance(G4FieldTrack& y_current,
		                 G4double     hstep,
		                 G4double     eps    )
//		     const G4double dydx[6],    // We could may add this ??

// Runge-Kutta driver with adaptive stepsize control. Integrate starting
// values at y_current over hstep x2 with accuracy eps. 
// On output ystart is replaced by values at the end of the integration 
// interval.  
// RightHandSide is the right-hand side of ODE system. 
// The source is similar to odeint routine from NRC p.721-722 .

{
  G4int nstp, i; 
  G4double x, hnext, hdid, h ;

  G4int  no_warnings=0;
#ifdef G4DEBUG_FIELD
  static G4int dbg=1;
  G4double ySubStepStart[G4FieldTrack::ncompSVEC];
  G4FieldTrack  yFldTrkStart(y_current);
#endif

  // G4double yscal[ncompSVEC];
  G4double y[G4FieldTrack::ncompSVEC], dydx[G4FieldTrack::ncompSVEC];
  G4double ystart[G4FieldTrack::ncompSVEC], yEnd[G4FieldTrack::ncompSVEC]; 
  G4double  x1, x2;
  G4bool succeeded = true, lastStepSucceeded;

  G4FieldTrack yStartFT(y_current);

  //  Assume that hstep > 0 

  // ystart = y_current.PosVelVec();
  y_current.DumpToArray( ystart );
  x1= y_current.GetCurveLength();
  x2= x1 + hstep;
 
  // //  Initial Step size "h" is half the interval
  // h = 0.5 * hstep;
  //  Initial Step size "h" is the full interval
  h = hstep;   
  x = x1;

  G4int  noFullIntegr=0, noSmallIntegr = 0 ;
  static G4int  noGoodSteps =0, noBadSteps = 0 ;  // Bad = chord > curve-len 

  for(i=0;i<nvar;i++) y[i] = ystart[i] ;

  G4bool  lastStep= false;
  nstp=1;
  // G4double  lastStepThreshold = G4std::min( eps * hstep, Hmin() ); 

  do{
     G4ThreeVector StartPos( y[0], y[1], y[2] );   
#    ifdef G4DEBUG_FIELD
       for(i=0;i<nvar;i++) ySubStepStart[i] = y[i] ;
       yFldTrkStart.LoadFromArray(y);
       yFldTrkStart.SetCurveLength(x);
#    endif

     pIntStepper->RightHandSide( y, dydx );

     if( x+h > x2 ) {
        h = x2 - x ;     // When stepsize overshoots, decrease it!
        if( h < eps * hstep) {
	    lastStep = true;   //  Ensure that this must be the last step
                               //   - since otherwise numerical (im)precision
                               //     could force lots of small last steps
	}
     }

     static G4int nStpPr=50;   // For debug printing of integrations with many steps

     // Perform the Integration
     //      

     if( h > Hmin() ){ 
        OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
        //--------------------------------------
        lastStepSucceeded= (hdid == h);   
#         ifdef  G4DEBUG_FIELD
  	  if(dbg>2) PrintStatus( ySubStepStart, x1, y, x, h,  nstp); // Only
#         endif
     }else{
        G4FieldTrack yFldTrk( G4ThreeVector(0,0,0), 
			      G4ThreeVector(0,0,0), 0., 0., 0., 0. );
        G4double dchord_step, dyerr, dyerr_len;  //  Must figure what to do with these
	yFldTrk.LoadFromArray(y); 
        yFldTrk.SetCurveLength( x );

        QuickAdvance( yFldTrk, dydx, h, dchord_step, dyerr_len ); 
        //-----------------------------------------------------

#         ifdef  G4DEBUG_FIELD
 	   // if(dbg>1) OneGoodStep(y,dydx,x,h,2*eps,hdid,hnext) ;
	   // if(dbg>1) PrintStatus( ystart, x1, y, x, h, -nstp);  
#         endif
        yFldTrk.DumpToArray(y);    
#         ifdef  G4DEBUG_FIELD
  	  if(dbg>1) PrintStatus( ySubStepStart, x1, y, x, h,  nstp);   // Only this
#         endif	
	dyerr = dyerr_len / h;    // was dyerr_len / hstep;
	hdid= h;
        x += hdid;
        // Compute suggested new step
	hnext= ComputeNewStepSize( dyerr/eps, h);
	// .. hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
	lastStepSucceeded= (dyerr<= eps);
     }

     if(lastStepSucceeded) noFullIntegr++ ; else noSmallIntegr++ ;
     G4ThreeVector EndPos( y[0], y[1], y[2] );

#ifdef  G4DEBUG_FIELD
     if(dbg && (nstp>nStpPr)) {
       G4cout << "hdid="  << G4std::setw(12) << hdid  << " "
	      << "hnext=" << G4std::setw(12) << hnext << " " << G4endl;
       PrintStatus( ystart, x1, y, x, h, (nstp==nStpPr) ? -nstp: nstp); 
     }
#endif

     // Check the endpoint
     G4double endPointDist= (EndPos-StartPos).mag(); 
     if( endPointDist >= hdid*(1.+perMillion) ){
        noBadSteps  ++;
        // Issue a warning only for gross differences -
        //   we understand how small difference occur.
        if( endPointDist >= hdid*(1.+perThousand) ){ 
#ifdef  G4DEBUG_FIELD
           if(dbg){
	      WarnEndPointTooFar ( endPointDist, hdid, eps, dbg ); 
	      G4cerr << "  Total steps:  bad" << noBadSteps << " good " << noGoodSteps << G4endl;
	      // G4cerr << "Mid:EndPtFar> ";
	      PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);  
                // Potentially add as arguments:  <dydx> - as Initial Force
           }
#endif
	   no_warnings++;
	}
     } else { // ie (!dbg)
        noGoodSteps ++;
     } 
// #endif

     // Check the proposed next stepsize
     if(fabs(hnext) <= Hmin())
     {
#ifdef  G4DEBUG_FIELD
        // If simply a very small interval is being integrated, do not warn
        if( (x < x2 * (1-eps) ) &&     //  The last step can be small: it's OK
            (fabs(hstep) > Hmin())     //   and if we are asked, it's OK
            //   && (hnext < hstep * PerThousand ) 
	  ){
             if(dbg>0){ 
		//  Issue WARNING
		WarnSmallStepSize( hnext, hstep, h, x-x1, nstp ); 
		// G4cerr << "Mid:SmallStep> ";
		PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);
	     }
	     no_warnings++;
	}
#endif
        // else 
	//   succeeded = false;  // Meaningful only if we break out of the loop.
	// 
        // lastStep = true;   //  Make this the last step ... Dubious now

        // Make sure that the next step is at least Hmin.
        h = Hmin();
     }else{
        h = hnext ;
     }

  }while ( ((nstp++)<=fMaxNoSteps)
          && (x < x2)           //  Have we reached the end ?
                                //   --> a better test might be x-x2 > an_epsilon
          && (!lastStep)
         );

  succeeded=  (x>=x2);  // If it was a "forced" last step

  for(i=0;i<nvar;i++)  yEnd[i] = y[i] ;

  // Put back the values.
  y_current.LoadFromArray( yEnd );
  y_current.SetCurveLength( x );

  if(nstp > fMaxNoSteps){
     no_warnings++;
     succeeded = false;
#ifdef  G4DEBUG_FIELD
        if(dbg){
	   WarnTooManyStep( x1, x2, x );  //  Issue WARNING
	   PrintStatus( yEnd, x1, y, x, hstep, -nstp);
        }
#endif
  }

#ifdef G4DEBUG_FIELD
  if( dbg && no_warnings ){
     G4cerr << " Exiting status: "
	   << " no-steps " << nstp
	   <<G4endl;
     PrintStatus( yEnd, x1, y, x, hstep, nstp);
  }
#endif

  return succeeded;

}  // end of AccurateAdvance ...........................

void
G4MagInt_Driver::WarnSmallStepSize( G4double hnext, G4double hstep, 
				    G4double h, G4double xDone,
				    G4int nstp)
{
  static G4int noWarningsIssued =0;
  const  G4int maxNoWarnings =  10;   // Number of verbose warnings
  if( (noWarningsIssued < maxNoWarnings) || fVerboseLevel > 10 ){
    G4cerr<< " Warning (G4MagIntegratorDriver::AccurateAdvance): The stepsize for the " 
	  << " next iteration=" << hnext << " is too small " 
	  << "- in Step number " << nstp << "." << G4endl;
    G4cerr << "     The minimum for the driver is " << Hmin()  << G4endl ;
    G4cerr << "     Requested integr. length was " << hstep << " ." << G4endl ;
    G4cerr << "     The size of this sub-step was " << h     << " ." << G4endl ;
    G4cerr << "     The integrations has already gone " << xDone << G4endl ;
  }else{
    G4cerr<< " G4MagInt_Driver: Too small 'next' step " << hnext 
	  << " step-no "  << nstp ;                     // << G4setw(4)
    G4cerr << "  this sub-step " << h     
	   << "  req_tot_len " << hstep 
	   << "  done " << xDone 
	   << "  min " << Hmin() 
	   << G4endl ;
  }
  noWarningsIssued++;
}

void
G4MagInt_Driver::WarnTooManyStep( G4double x1start, 
				  G4double x2end, 
				  G4double xCurrent)
{
    G4cerr << " Warning (G4MagIntegratorDriver): The number of steps " 
         << "used in the Integration driver (Runge-Kutta) is too many.  "
	 << G4endl ;
    G4cerr << "Integration of the interval was not completed - only a " 
         << (xCurrent-x1start)*100/(x2end-x1start)
	   <<" % fraction of it was Done." << G4endl;
}

void
G4MagInt_Driver::WarnEndPointTooFar (G4double endPointDist, 
				     G4double   h , 
				     G4double  eps,
				     G4int     dbg)
{
	static G4double maxRelError= 0.0, maxRelError_last_printed=0.0;
	G4bool isNewMax, prNewMax;

        isNewMax = endPointDist > (1.0 + maxRelError) * h;
        prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
	if( isNewMax )
	   maxRelError= endPointDist / h - 1.0; 
        if( prNewMax )
 	   maxRelError_last_printed = maxRelError;

        if(    dbg 
	    && (h > kCarTolerance) 
	    && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+eps) ) ) 
          ){ 
           static G4int noWarnings = 0;
           if( (noWarnings ++ < 10) || (dbg>2) ){
 	      G4cerr << " Warning (G4MagIntegratorDriver): "
		     << " The integration produced an endpoint which " << G4endl
		     << "   is further from the startpoint than the curve length." << G4endl; 
          
	      G4cerr << "   Distance of endpoints = " << endPointDist
		     << "  curve length = " <<  h
		     << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
		     << "  relative = " << (h-endPointDist) / h 
		     << "  epsilon =  " << eps 
		     << G4endl;
	   }else{
 	      G4cerr << " Warning:" 
		     << "  dist_e= " << endPointDist
		     << "  h_step = " <<  h
		     << "  Diff (hs-ed)= " << (h - endPointDist)
		     << "  rel = " << (h-endPointDist) / h 
		     << "  eps = " << eps
		     << " (from G4MagInt_Driver)" << G4endl;
	   }
	}
}
// ---------------------------------------------------------

void
G4MagInt_Driver::OneGoodStep(      G4double y[],        // InOut
			     const G4double dydx[],
				   G4double& x,         // InOut
			           G4double htry,
			           G4double eps_rel_max,
				   G4double& hdid,      // Out
				   G4double& hnext )    // Out

// Driver for one Runge-Kutta Step with monitoring of local truncation error
// to ensure accuracy and adjust stepsize. Input are dependent variable
// array y[0,...,5] and its derivative dydx[0,...,5] at the
// starting value of the independent variable x . Also input are stepsize
// to be attempted htry, and the required accuracy eps. On output y and x
// are replaced by their new values, hdid is the stepsize that was actually
// accomplished, and hnext is the estimated next stepsize. 
// This is similar to the function rkqs from the book:
// Numerical Recipes in C: The Art of Scientific Computing (NRC), Second
// Edition, by William H. Press, Saul A. Teukolsky, William T.
// Vetterling, and Brian P. Flannery (Cambridge University Press 1992),
// 16.2 Adaptive StepSize Control for Runge-Kutta, p. 719

{
      G4double errpos_sq, errvel_sq, errmax_sq;
      G4double errmax, h, htemp, xnew ;
      G4int i;

      G4double yerr[G4FieldTrack::ncompSVEC], ytemp[G4FieldTrack::ncompSVEC];

      h = htry ; // Set stepsize to the initial trial value

      // G4double inv_epspos_sq= 1.0 / eps * eps; 

      for (;;)
      {
	  pIntStepper-> Stepper(y,dydx,h,ytemp,yerr); 
          G4double eps_pos = eps_rel_max * G4std::max(h, Hmin()); 
	  // Evaluate accuracy
	  //
	  errpos_sq =  sqr(yerr[0]) + sqr(yerr[1]) + sqr(yerr[2]) ;
	  errpos_sq /= eps_pos*eps_pos; // Scale relative to required tolerance

          // Accuracy for velocity
          errvel_sq =  (sqr(yerr[3]) + sqr(yerr[4]) + sqr(yerr[5]) )
                     / (sqr(y[3]) + sqr(y[4]) + sqr(y[5]) );
          errvel_sq /= eps_rel_max*eps_rel_max; 

          errmax_sq = G4std::max( errpos_sq, errvel_sq ); // Square of maximum error
          errmax = sqrt( errmax_sq );
	  if(errmax_sq <= 1.0 ) break ; // Step succeeded. 

	  // Step failed; compute the size of retrial Step.
	  htemp = GetSafety()*h*pow(errmax,GetPshrnk()) ;

	  if(htemp >= 0.1*h) h = htemp ;  // Truncation error too large,
	  else h = 0.1*h ;                // reduce stepsize, but no more
					  // than a factor of 10
	  xnew = x + h ;
	  if(xnew == x) {
	     G4cerr<<"G4MagIntegratorDriver::OneGoodStep: Stepsize underflow in Stepper "<<G4endl ;
	     G4cerr<<"  Step's start x=" << x << " and end x= " << xnew 
		   << " are equal !! " << G4endl
		   <<"  Due to step-size= " << h 
                   << " . Note that input step was " << htry << G4endl;
	     break;
	  }
      }

      // Compute size of next Step
      if(errmax > errcon) hnext = GetSafety()*h*pow(errmax,GetPgrow()) ;
      else hnext = max_stepping_increase*h ;
                     // No more than a factor of 5 increase

      x += (hdid = h) ;

      for(i=0;i<nvar;i++) y[i] = ytemp[i] ;

      // delete[] ytemp ;
      // delete[] yerr  ;
      return ;

}   // end of  OneGoodStep .............................


//----------------------------------------------------------------------
// QuickAdvance just tries one Step - it does not ensure accuracy
//
G4bool  G4MagInt_Driver::QuickAdvance(       
			    G4FieldTrack& y_posvel,         // INOUT
		            const G4double     dydx[],  
		                  G4double     hstep,       // In
			          G4double&    dchord_step,
			          G4double&    dyerr )  
{
    G4double yerr_vec[G4FieldTrack::ncompSVEC], yarrin[G4FieldTrack::ncompSVEC], yarrout[G4FieldTrack::ncompSVEC]; 
    G4double s_start;
    G4double dyerr_len, dyerr_vel, vel_mag;

    // Move data into array
    y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
    s_start = y_posvel.GetCurveLength();

    // Do an Integration Step
    pIntStepper-> Stepper(yarrin, dydx, hstep, yarrout, yerr_vec) ; 

    // Estimate curve-chord distance
    dchord_step= pIntStepper-> DistChord();

    // Put back the values.
    y_posvel.LoadFromArray( yarrout );   //  yarrout ==> y_posvel
    y_posvel.SetCurveLength( s_start + hstep );

    // A single measure of the error   
    //      TO-DO :  account for  tangent vector,  energy,  spin, ... ? 
    dyerr_len= sqrt( sqr(yerr_vec[0])+sqr(yerr_vec[1])+sqr(yerr_vec[2]));
    dyerr_vel= sqrt( sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));
    vel_mag  = sqrt( sqr(yarrout[3])+sqr(yarrout[4])+sqr(yarrout[5]) );

    if( (dyerr_len / hstep) > (dyerr_vel / vel_mag) ) {
       dyerr = dyerr_len;
    }else{
       // Scale it to the position - for now
       dyerr = (dyerr_vel / vel_mag) * hstep;
    }
#ifdef RETURN_A_NEW_STEP_LENGTH
    // The following step cannot be done here because "eps" is not known.
    dyerr_len /= eps;

    // Look at the velocity deviation ?
    //  sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));

    // Look at the change in the velocity (squared maybe ..)
    G4double veloc_square = y_posvel.GetVelocity().mag2();

    // Set suggested new step
    hstep= ComputeNewStepSize( dyerr_len, hstep);
#endif

    return true;
}

#ifdef QUICK_ADV_TWO
G4bool  G4MagInt_Driver::QuickAdvance(       
				  G4double     yarrin[],        // IN
		            const G4double     dydx[],  
		                  G4double     hstep,       // In
				  G4double     yarrout[],
			          G4double&    dchord_step,
			          G4double&    dyerr )      // in length
{
   G4Exception("Not implemented in current version");

   dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
   yarrout[0]= yarrin[0];
}
#endif 

// --------------------------------------------------------------------------
//  This method computes new step sizes - but does not limit changes to
//   within  certain factors
// 

G4double 
G4MagInt_Driver::ComputeNewStepSize( 
                          G4double  errMaxNorm,    // max error  (normalised)
			  G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if(errMaxNorm > 1.0 ) {

    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*pow(errMaxNorm,GetPshrnk()) ;
  }else if(errMaxNorm > 0.0 ){
    // Compute size of next Step for a successful step
    hnew = GetSafety()*hstepCurrent*pow(errMaxNorm,GetPgrow()) ;
  }else {
    // if error estimate is zero (possible) or negative (dubious)
    hnew = max_stepping_increase * hstepCurrent; 
  }

  return hnew;
}

// -----------------------------------------------------------------------------
//  This method computes new step sizes limiting changes within certain factors
// 
//   It shares its logic with AccurateAdvance.
//    They are kept separate currently for optimisation.

G4double 
G4MagInt_Driver::ComputeNewStepSize_WithinLimits( 
                          G4double  errMaxNorm,    // max error  (normalised)
			  G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if(errMaxNorm > 1.0 ) {

    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*pow(errMaxNorm,GetPshrnk()) ;
  
    if(hnew < max_stepping_decrease*hstepCurrent) 
         hnew = max_stepping_decrease*hstepCurrent ;
                         // reduce stepsize, but no more
                         // than this factor (value= 1/10)
  }else{
    // Compute size of next Step for a successful step
    if(errMaxNorm > errcon) hnew = GetSafety()*hstepCurrent*pow(errMaxNorm,GetPgrow()) ;
    else                    hnew = max_stepping_increase * hstepCurrent ;
      // No more than a factor of 5 increase
  }

  return hnew;
}



void G4MagInt_Driver::PrintStatus( const G4double*   StartArr,  
				   G4double          xstart,
				   const G4double*   CurrentArr, 
				   G4double          xcurrent,
				   G4double          requestStep, 
				   G4int             subStepNo)
  // Potentially add as arguments:  
  //                                 <dydx>           - as Initial Force
  //                                 stepTaken(hdid)  - last step taken
  //                                 nextStep (hnext) - proposal for size
{
   G4FieldTrack  StartFT(G4ThreeVector(0,0,0), G4ThreeVector(0,0,0), 0., 0., 0., 0. );
   G4FieldTrack  CurrentFT (StartFT);

   StartFT.LoadFromArray( StartArr); 
   StartFT.SetCurveLength( xstart);
   CurrentFT.LoadFromArray( CurrentArr); 
   CurrentFT.SetCurveLength( xcurrent );

   PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
}



void G4MagInt_Driver::PrintStatus(
                  const G4FieldTrack&  StartFT,
		  const G4FieldTrack&  CurrentFT, 
                  G4double             requestStep, 
                  // G4double             safety,
                  G4int                subStepNo)
{
    G4int verboseLevel= fVerboseLevel;
    static G4int noPrecision= 5;
    G4int oldPrec= G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);

    const G4ThreeVector StartPosition=      StartFT.GetPosition();
    const G4ThreeVector StartUnitVelocity=  StartFT.GetMomentumDir();
    const G4ThreeVector CurrentPosition=    CurrentFT.GetPosition();
    const G4ThreeVector CurrentUnitVelocity=    CurrentFT.GetMomentumDir();

    G4double  DotStartCurrentVeloc= StartUnitVelocity.dot(CurrentUnitVelocity);

    G4double step_len= CurrentFT.GetCurveLength() 
                         - StartFT.GetCurveLength();
    G4double subStepSize = step_len;
     
    if( (subStepNo <= 0) && (verboseLevel <= 3) )
    {
       subStepNo = - subStepNo;        // To allow printing banner

       G4cout << G4std::setw( 6)  << " " 
	      << G4std::setw( 25) << " G4MagInt_Driver: Current Position  and  Direction" << " "
	      << G4endl; 
       G4cout << G4std::setw( 5) << "Step#" << " "
	      << G4std::setw( 7) << "s-curve" << " "
	      << G4std::setw( 9) << "X(mm)" << " "
	      << G4std::setw( 9) << "Y(mm)" << " "  
	      << G4std::setw( 9) << "Z(mm)" << " "
	      << G4std::setw( 8) << " N_x " << " "
	      << G4std::setw( 8) << " N_y " << " "
	      << G4std::setw( 8) << " N_z " << " "
	      << G4std::setw( 7) << " N^2-1 " << " "
	      << G4std::setw(10) << " N(0).N " << " "
	      << G4std::setw( 7) << "KinEner " << " "
	      << G4std::setw(12) << "Track-l" << " "   // Add the Sub-step ??
	      << G4std::setw(12) << "Step-len" << " " 
	      << G4std::setw(12) << "Step-len" << " " 
	      << G4std::setw( 9) << "ReqStep" << " "  
	      << G4endl;

        PrintStat_Aux( StartFT,  requestStep, 0., 
		       0,        0.0,         1.0);
        //*************
    }

    if( verboseLevel <= 3 )
    {
       G4cout.precision(noPrecision);
       PrintStat_Aux( CurrentFT, requestStep, step_len, 
		      subStepNo, subStepSize, DotStartCurrentVeloc );
       //*************
    }

    else // if( verboseLevel > 3 )
    {
       //  Multi-line output
       
       // G4cout << "Current  Position is " << CurrentPosition << G4endl 
       //    << " and UnitVelocity is " << CurrentUnitVelocity << G4endl;
       // G4cout << "Step taken was " << step_len  
       //    << " out of PhysicalStep= " <<  requestStep << G4endl;
       // G4cout << "Final safety is: " << safety << G4endl;

       // G4cout << "Chord length = " << (CurrentPosition-StartPosition).mag() << G4endl;
       // G4cout << G4endl; 
    }
    G4cout.precision(oldPrec);
}

void G4MagInt_Driver::PrintStat_Aux(
                  const G4FieldTrack&  aFieldTrack,
                  G4double             requestStep, 
		  G4double             step_len,
                  G4int                subStepNo,
		  G4double             subStepSize,
		  G4double             dotVeloc_StartCurr)
{
    const G4ThreeVector Position=      aFieldTrack.GetPosition();
    const G4ThreeVector UnitVelocity=  aFieldTrack.GetMomentumDir();
 
    if( subStepNo >= 0)
       G4cout << G4std::setw( 5) << subStepNo << " ";
    else
       G4cout << G4std::setw( 5) << "Start" << " ";
    G4double curveLen= aFieldTrack.GetCurveLength();
    G4cout << G4std::setw( 7) << curveLen;
    G4cout << G4std::setw( 9) << Position.x() << " "
	   << G4std::setw( 9) << Position.y() << " "
	   << G4std::setw( 9) << Position.z() << " "
	   << G4std::setw( 8) << UnitVelocity.x() << " "
	   << G4std::setw( 8) << UnitVelocity.y() << " "
	   << G4std::setw( 8) << UnitVelocity.z() << " ";
    G4int oldprec= G4cout.precision(3);
    G4cout << G4std::setw( 7) << UnitVelocity.mag2()-1.0 << " ";
    G4cout.precision(6);
    G4cout << G4std::setw(10) << dotVeloc_StartCurr << " ";
    G4cout.precision(oldprec);
    G4cout << G4std::setw( 7) << aFieldTrack.GetKineticEnergy();
    G4cout << G4std::setw(12) << step_len << " ";

    static G4double oldCurveLength= 0.0;
    static G4double oldSubStepLength= 0.0;
    static int oldSubStepNo= -1;

    G4double subStep_len=0.0;
    if( curveLen > oldCurveLength )
       subStep_len= curveLen - oldCurveLength;
    else if (subStepNo == oldSubStepNo)
       subStep_len= oldSubStepLength;
    //     else  subStepLen_NotAvail;
    oldCurveLength= curveLen;
    oldSubStepLength= subStep_len;

    G4cout << G4std::setw(12) << subStep_len << " "; 
    G4cout << G4std::setw(12) << subStepSize << " "; 
    if( requestStep != -1.0 ) 
       G4cout << G4std::setw( 9) << requestStep << " ";
    else
       G4cout << G4std::setw( 9) << " InitialStep " << " "; 
    // G4cout << G4std::setw(12) << safety << " ";
    G4cout << G4endl;
}
