// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagIntegratorDriver.cc,v 1.11 2000-05-05 18:15:49 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Implementation for class G4MagInt_Driver
//   Tracking in space dependent magnetic field
//
// History of major changes:
//  7 Oct 96  V. Grichine       First version
// 28 Jan 98  W. Wander:        Added ability for low order integrators
// 30 Jan 98  J. Apostolakis:   Made method parameters into instance variables
// 27 Jul 99  J. Apostolakis:   Ensured that AccurateAdvance does not loop 
//                                due to very small eps & step size (precision)
#include <math.h>
#include "G4ios.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4FieldTrack.hh"
#include "geomdefs.hh"          
                             //  for kCarTolerance

#define G4DEBUG 1

//  Stepsize can increase by no more than 5.0
//           and decrease by no more than 1/10. = 0.1
//
const G4double G4MagInt_Driver::max_stepping_increase = 5.0;
const G4double G4MagInt_Driver::max_stepping_decrease = 0.1;

//  The (default) maximum number of steps is Base divided by the order of Stepper
const G4int  G4MagInt_Driver::fMaxStepBase = 5000;  

G4bool
G4MagInt_Driver::AccurateAdvance( 
			   G4FieldTrack& y_current,
		     const G4double     hstep,
		     const G4double     eps
		     )
//		     const G4double dydx[6],    // We could may add this ??

// Runge-Kutta driver with adaptive stepsize control. Integrate starting
// values at y_current over hstep x2 with accuracy eps. 
// On output ystart is replaced by values at the end of the integration 
// interval.  
// RightHandSide is the right-hand side of ODE system. 
// The source is similar to odeint routine from NRC p.721-722 .

// OLD:
// The value h1 should be set as a guessed first stepsize, Hmin is the
// minimum allowed stepsize. On output nOK and nBAD are the numbers of
// good and bad (but retried and fixed) steps taken.
{
  // static const G4int maxstp = 5000;

  G4int nstp, i; 
  static G4int dbg=1;
  G4double x, hnext, hdid, h ;

  // G4double yscal[ncompSVEC];
  G4double y[G4FieldTrack::ncompSVEC], dydx[G4FieldTrack::ncompSVEC];
  G4double ystart[G4FieldTrack::ncompSVEC]; 
  G4double  x1, x2;
  G4bool succeeded = true, lastStepSucceeded;

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
  G4double  lastStepThreshold = G4std::min( eps * hstep, Hmin() ); 

  do{
     G4ThreeVector StartPos( y[0], y[1], y[2] );   

     pIntStepper->RightHandSide( y, dydx );

     if( x+h > x2 ) {
       h = x2 - x ;     // When stepsize overshoots, decrease it!
     }
     if( h < eps * hstep) {
       lastStep = true;   //  Ensure that this must be the last step
                          //   because otherwise numerical (im)precision
                          //   could otherwise force lots of small last steps.
     }

     // Perform the Integration
     // 
     OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
     //--------------------------------------

     lastStepSucceeded= (hdid == h);   
#ifdef G4DEBUG
     if(lastStepSucceeded) noFullIntegr++ ; else noSmallIntegr++ ;
     G4ThreeVector EndPos( y[0], y[1], y[2] );

     // Check the endpoint
     G4double endPointDist= (EndPos-StartPos).mag(); 
     if( endPointDist >= h*(1.+perMillion) ){
        // Issue a warning only for gross differences -
        //   we understand how small difference occur.
        if( endPointDist >= h*(1.+perThousand) ){ 
           WarnEndPointTooFar ( endPointDist, h, eps, dbg ); 
	}
        noBadSteps  ++;
     } else { // ie (!dbg)
        noGoodSteps ++;
     } 
     
#endif

     // Check the proposed next stepsize
     if(fabs(hnext) <= Hmin())
     {
        // If simply a very small interval is being integrated, do not warn
        if( (x < x2 * (1-eps) ) &&     //  The last step can be small: it's OK
            (fabs(hstep) > Hmin())     //   and if we are asked, it's OK
            //   && (hnext < hstep * PerThousand ) 
          )
	  //  Issue WARNING
	  WarnSmallStepSize( hnext, hstep, h, x-x1, nstp ); 
        else 
	  succeeded = false;  // Meaningful only if we break out of the loop.

        lastStep = true;   // ensure that this was the last step
     }

     h = hnext ;
  }while (((nstp++)<=fMaxNoSteps) &&
          (x < x2)           //  Have we reached the end ?
                             //   --> a better test might be x-x2 > an_epsilon
          && (!lastStep)
         );

  succeeded=  (x>=x2);  // If it was a "forced" last step

  for(i=0;i<nvar;i++)  ystart[i] = y[i] ;

  if(nstp > fMaxNoSteps){
     WarnTooManyStep( x1, x2, x );  //  Issue WARNING
     succeeded = false;
  }

  // Put back the values.
  y_current.LoadFromArray( ystart );
  y_current.SetCurveLength( x );

  return succeeded;

}  // end of AccurateAdvance ...........................

void
G4MagInt_Driver::WarnSmallStepSize( G4double hnext, G4double hstep, 
				    G4double h, G4double xDone,
				    G4int nstp)
{
  static G4int noWarningsIssued =0;
  const  G4int maxNoWarnings = 100; 
  if( noWarningsIssued < maxNoWarnings ){
    G4cerr<< " Warning (G4MagIntegratorDriver): The stepsize for the " 
	  << " next iteration=" << hnext << " is too small - in Step number "
	  << nstp << "." << G4endl;
    G4cerr << "     Requested step size was " << hstep << " ." << G4endl ;
    G4cerr << "     Previous  step size was " << h     << " ." << G4endl ;
    G4cerr << "     The minimum for the driver is " << Hmin()  << G4endl ;
    G4cerr << "     The integrations has already gone " << xDone << G4endl ;
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
	static G4double maxRelError= 0.0;
	G4bool isNewMax, prNewMax;

        isNewMax = endPointDist > (1.0 + maxRelError) * h;
        prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
	if( isNewMax )
	   maxRelError= endPointDist / h - 1.0; 

        if(    dbg 
	    && (h > kCarTolerance) 
	    && ( prNewMax || (endPointDist >= h*(1.+eps) ) ) 
          ){ 
           static G4int noWarnings = 0;
           if( (noWarnings ++ < 10) || (dbg>1) ){
 	      G4cerr << " Warning (G4MagIntegratorDriver): "
		     << " The integration produced an endpoint which " << G4endl
		     << "   is further from the startpoint than the curve length." << G4endl; 
          
	      G4cerr << "   Distance of endpoints = " << endPointDist
		     << "  curve length = " <<  h
		     << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
		     << "  relative = " << (h-endPointDist) / h 
		     << G4endl;
	   }else{
	      G4cerr << "  EndpointDist = " << endPointDist
		     << "  curve length = " <<  h
		     << "  Diff (cl-ed)= " << (h - endPointDist)
		     << "  rel = " << (h-endPointDist) / h 
		     << " (from G4MagIntegratorDriver)" << G4endl;
	   }
	}
}
// ---------------------------------------------------------

void
G4MagInt_Driver::OneGoodStep(      G4double y[],
			     const G4double dydx[],
				   G4double& x,
			     const G4double htry,
			     const G4double eps_rel_max,
				   G4double& hdid,
				   G4double& hnext )

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
			    G4FieldTrack& y_posvel,   // INOUT
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
  }else{
    // Compute size of next Step for a successful step
    hnew = GetSafety()*hstepCurrent*pow(errMaxNorm,GetPgrow()) ;
  }

  return hnew;
}

// --------------------------------------------------------------------------
//  This method computes new step sizes - but does not limit changes to
//   within  certain factors
// 
//   It shares its logic with AccurateAdvance, so they should eventually 
//     be merged  ??

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
