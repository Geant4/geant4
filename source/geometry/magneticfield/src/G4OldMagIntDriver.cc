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
// G4OldMagIntDriver -- same behaviour as old G4MagInt_Driver
//
// V.Grichine, 07.10.1996 - Created
// J.Apostolakis, 08.11.2001 - Respect minimum step in AccurateAdvance
// --------------------------------------------------------------------

#include <iomanip>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"
#include "G4OldMagIntDriver.hh"
#include "G4FieldTrack.hh"

#ifdef   G4DEBUG_FIELD
#include "G4DriverReporter.hh"
#endif

// ---------------------------------------------------------

//  Constructor
//
G4OldMagIntDriver::G4OldMagIntDriver( G4double                hminimum, 
                                  G4MagIntegratorStepper* pStepper,
                                  G4int                   numComponents,
                                  G4int                   statisticsVerbose)
  : fNoIntegrationVariables(numComponents), 
    fNoVars( std::max( fNoIntegrationVariables, fMinNoVars )),
    fStatisticsVerboseLevel(statisticsVerbose)
{  
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12

  RenewStepperAndAdjust( pStepper );
  fMinimumStep = hminimum;

  fMaxNoSteps = fMaxStepBase / pIntStepper->IntegratorOrder();
#ifdef G4DEBUG_FIELD
  fVerboseLevel=2;
#endif

  if( (fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1) )
  {
    G4cout << "MagIntDriver version: Accur-Adv: "
           << "invE_nS, QuickAdv-2sqrt with Statistics "
#ifdef G4FLD_STATS
           << " enabled "
#else
           << " disabled "
#endif
           << G4endl;
  }
}

// ---------------------------------------------------------

//  Destructor
//
G4OldMagIntDriver::~G4OldMagIntDriver()
{ 
  if( fStatisticsVerboseLevel > 1 )
  {
    PrintStatisticsReport();
  }
}

// ---------------------------------------------------------

G4bool
G4OldMagIntDriver::AccurateAdvance(G4FieldTrack& y_current,
                                 G4double      hstep,
                                 G4double      eps,
                                 G4double      hinitial )
{
  // Runge-Kutta driver with adaptive stepsize control. Integrate starting
  // values at y_current over hstep x2 with accuracy eps. 
  // On output ystart is replaced by values at the end of the integration 
  // interval. RightHandSide is the right-hand side of ODE system. 
  // The source is similar to odeint routine from NRC p.721-722 .

  G4int nstp, i;
  G4double x, hnext, hdid, h;

#ifdef G4DEBUG_FIELD
  G4int no_warnings = 0;
  static G4int dbg = 1;
  G4double ySubStepStart[G4FieldTrack::ncompSVEC];
  G4FieldTrack  yFldTrkStart(y_current);
#endif

  G4double y[G4FieldTrack::ncompSVEC], dydx[G4FieldTrack::ncompSVEC];
  G4double ystart[G4FieldTrack::ncompSVEC], yEnd[G4FieldTrack::ncompSVEC]; 
  G4double  x1, x2;
  G4bool succeeded = true;

  G4double startCurveLength;

  const G4int nvar = fNoVars;

  G4FieldTrack yStartFT(y_current);

  //  Ensure that hstep > 0
  //
  if( hstep <= 0.0 )
  { 
    if( hstep == 0.0 )
    {
      std::ostringstream message;
      message << "Proposed step is zero; hstep = " << hstep << " !";
      G4Exception("G4OldMagIntDriver::AccurateAdvance()", 
                  "GeomField1001", JustWarning, message);
      return succeeded; 
    }
    else
    { 
      std::ostringstream message;
      message << "Invalid run condition." << G4endl
              << "Proposed step is negative; hstep = " << hstep << "." << G4endl
              << "Requested step cannot be negative! Aborting event.";
      G4Exception("G4OldMagIntDriver::AccurateAdvance()", 
                  "GeomField0003", EventMustBeAborted, message);
      return false;
    }
  }

  y_current.DumpToArray( ystart );

  startCurveLength= y_current.GetCurveLength();
  x1= startCurveLength; 
  x2= x1 + hstep;

  if ( (hinitial > 0.0) && (hinitial < hstep)
    && (hinitial > perMillion * hstep) )
  {
     h = hinitial;
  }
  else  //  Initial Step size "h" defaults to the full interval
  {
     h = hstep;
  }

  x = x1;

  for ( i=0; i<nvar; ++i)  { y[i] = ystart[i]; }

#ifdef G4DEBUG_FIELD  
  // G4cout << "IDriver: hstep = " << hstep << " hinitial= " << hinitial << " h = " << h << G4endl;
  G4cout << "IDriver::AccurAdv called. " 
         << " Input: hstep = " << hstep << " hinitial= " << hinitial << " , current: h = " << h << G4endl;        
#endif
  
  G4bool lastStep= false;
  nstp = 1;

  do
  {
    G4ThreeVector StartPos( y[0], y[1], y[2] );

#ifdef G4DEBUG_FIELD
    G4double xSubStepStart= x; 
    for (i=0; i<nvar; ++i)  { ySubStepStart[i] = y[i]; }
    yFldTrkStart.LoadFromArray(y, fNoIntegrationVariables);
    yFldTrkStart.SetCurveLength(x);
    if(dbg) // Debug
       G4cout << "----- Iteration = " << nstp-1 << G4endl;
#endif

    pIntStepper->RightHandSide( y, dydx );
    ++fNoTotalSteps;

    // Perform the Integration
    //      
    if( h > fMinimumStep )
    { 
      OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
      //--------------------------------------

#ifdef G4DEBUG_FIELD
      if (dbg) // (dbg>2)
      {
        G4cout << "IntegrationDriver -- after OneGoodStep / requesting step = " << h << G4endl;
        // PrintStatus( ySubStepStart, xSubStepStart, y, x, h,  nstp); // Only        
        G4DriverReporter::PrintStatus( ySubStepStart, xSubStepStart, y, x, h,  nstp, nvar);
      }
#endif
    }
    else
    {
      G4FieldTrack yFldTrk( G4ThreeVector(0,0,0), 
                            G4ThreeVector(0,0,0), 0., 0., 0., 0. );
      G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
      yFldTrk.LoadFromArray(y, fNoIntegrationVariables); 
      yFldTrk.SetCurveLength( x );

      QuickAdvance( yFldTrk, dydx, h, dchord_step, dyerr_len );
      //-----------------------------------------------------

      yFldTrk.DumpToArray(y);    

#ifdef G4FLD_STATS
      ++fNoSmallSteps; 
      if ( dyerr_len > fDyerr_max )  { fDyerr_max = dyerr_len; }
      fDyerrPos_smTot += dyerr_len;
      fSumH_sm += h;  // Length total for 'small' steps
      if (nstp<=1)  { ++fNoInitialSmallSteps; }
#endif
#ifdef G4DEBUG_FIELD
      if (dbg>1)
      {
        if(fNoSmallSteps<2) { PrintStatus(ySubStepStart, x1, y, x, h, -nstp); }
        G4cout << "Another sub-min step, no " << fNoSmallSteps 
               << " of " << fNoTotalSteps << " this time " << nstp << G4endl; 
        PrintStatus( ySubStepStart, x1, y, x, h,  nstp);   // Only this
        G4cout << " dyerr= " << dyerr_len << " relative = " << dyerr_len / h 
               << " epsilon= " << eps << " hstep= " << hstep 
               << " h= " << h << " hmin= " << fMinimumStep << G4endl;
      }
#endif        
      if( h == 0.0 )
      { 
        G4Exception("G4OldMagIntDriver::AccurateAdvance()",
                    "GeomField0003", FatalException,
                    "Integration Step became Zero!"); 
      }
      dyerr = dyerr_len / h;
      hdid = h;
      x += hdid;

      // Compute suggested new step
      hnext = ComputeNewStepSize( dyerr/eps, h);

      // .. hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
    }

    G4ThreeVector EndPos( y[0], y[1], y[2] );

#if (G4DEBUG_FIELD>1)
    static G4int nStpPr = 50;   // For debug printing of long integrations
    if( (dbg>0) && (dbg<=2) && (nstp>nStpPr))
    {
      if( nstp==nStpPr )  { G4cout << "***** Many steps ****" << G4endl; }
      G4cout << "MagIntDrv: " ; 
      G4cout << "hdid="  << std::setw(12) << hdid  << " "
             << "hnext=" << std::setw(12) << hnext << " " 
             << "hstep=" << std::setw(12) << hstep << " (requested) " 
             << G4endl;
      PrintStatus( ystart, x1, y, x, h, (nstp==nStpPr) ? -nstp: nstp); 
    }
#endif

    // Check the endpoint
    G4double endPointDist= (EndPos-StartPos).mag(); 
    if ( endPointDist >= hdid*(1.+perMillion) )
    {
      ++fNoBadSteps;

      // Issue a warning only for gross differences -
      // we understand how small difference occur.
      if ( endPointDist >= hdid*(1.+perThousand) )
      { 
#ifdef G4DEBUG_FIELD
        if (dbg)
        {
          WarnEndPointTooFar ( endPointDist, hdid, eps, dbg ); 
          G4cerr << "  Total steps:  bad " << fNoBadSteps
                 << " current h= " << hdid << G4endl;
          PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);  
        }
        ++no_warnings;
#endif
      }
    }

    //  Avoid numerous small last steps
    if( (h < eps * hstep) || (h < fSmallestFraction * startCurveLength) )
    {
      // No more integration -- the next step will not happen
      lastStep = true;  
    }
    else
    {
      // Check the proposed next stepsize
      if(std::fabs(hnext) <= Hmin())
      {
#ifdef  G4DEBUG_FIELD
        // If simply a very small interval is being integrated, do not warn
        if( (x < x2 * (1-eps) ) &&        // The last step can be small: OK
            (std::fabs(hstep) > Hmin()) ) // and if we are asked, it's OK
        {
          if(dbg>0)
          {
            WarnSmallStepSize( hnext, hstep, h, x-x1, nstp );  
            PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);
          }
          ++no_warnings;
        }
#endif
        // Make sure that the next step is at least Hmin.
        h = Hmin();
      }
      else
      {
        h = hnext;
      }

      //  Ensure that the next step does not overshoot
      if ( x+h > x2 )
      {                // When stepsize overshoots, decrease it!
        h = x2 - x ;   // Must cope with difficult rounding-error
      }                // issues if hstep << x2

      if ( h == 0.0 )
      {
        // Cannot progress - accept this as last step - by default
        lastStep = true;
#ifdef G4DEBUG_FIELD
        if (dbg>2)
        {
          int prec= G4cout.precision(12); 
          G4cout << "Warning: G4MagIntegratorDriver::AccurateAdvance"
                 << G4endl
                 << "  Integration step 'h' became "
                 << h << " due to roundoff. " << G4endl
                 << " Calculated as difference of x2= "<< x2 << " and x=" << x
                 << "  Forcing termination of advance." << G4endl;
          G4cout.precision(prec);
        }          
#endif
      }
    }
  } while ( ((nstp++)<=fMaxNoSteps) && (x < x2) && (!lastStep) );
  // Loop checking, 07.10.2016, J. Apostolakis

     // Have we reached the end ?
     // --> a better test might be x-x2 > an_epsilon

  succeeded = (x>=x2);  // If it was a "forced" last step

  for (i=0; i<nvar; ++i)  { yEnd[i] = y[i]; }

  // Put back the values.
  y_current.LoadFromArray( yEnd, fNoIntegrationVariables );
  y_current.SetCurveLength( x );

  if(nstp > fMaxNoSteps)
  {
    succeeded = false;
#ifdef G4DEBUG_FIELD
    ++no_warnings;
    if (dbg)
    {
      WarnTooManyStep( x1, x2, x );  //  Issue WARNING
      PrintStatus( yEnd, x1, y, x, hstep, -nstp);
    }
#endif
  }

#ifdef G4DEBUG_FIELD
  if( dbg && no_warnings )
  {
    G4cerr << "G4MagIntegratorDriver exit status: no-steps " << nstp << G4endl;
    PrintStatus( yEnd, x1, y, x, hstep, nstp);
  }
#endif

  return succeeded;
}  // end of AccurateAdvance ...........................

// ---------------------------------------------------------

void
G4OldMagIntDriver::WarnSmallStepSize( G4double hnext, G4double hstep, 
                                    G4double h, G4double xDone,
                                    G4int nstp)
{
  static G4ThreadLocal G4int noWarningsIssued = 0;
  const  G4int maxNoWarnings = 10;   // Number of verbose warnings
  std::ostringstream message;
  if( (noWarningsIssued < maxNoWarnings) || fVerboseLevel > 10 )
  {
    message << "The stepsize for the next iteration, " << hnext
            << ", is too small - in Step number " << nstp << "." << G4endl
            << "The minimum for the driver is " << Hmin()  << G4endl
            << "Requested integr. length was " << hstep << " ." << G4endl
            << "The size of this sub-step was " << h     << " ." << G4endl
            << "The integrations has already gone " << xDone;
  }
  else
  {
    message << "Too small 'next' step " << hnext
            << ", step-no: " << nstp << G4endl
            << ", this sub-step: " << h     
            << ",  req_tot_len: " << hstep 
            << ", done: " << xDone << ", min: " << Hmin();
  }
  G4Exception("G4OldMagIntDriver::WarnSmallStepSize()", "GeomField1001",
              JustWarning, message);
  ++noWarningsIssued;
}

// ---------------------------------------------------------

void
G4OldMagIntDriver::WarnTooManyStep( G4double x1start, 
                                  G4double x2end, 
                                  G4double xCurrent )
{
   std::ostringstream message;
   message << "The number of steps used in the Integration driver"
           << " (Runge-Kutta) is too many." << G4endl
           << "Integration of the interval was not completed !" << G4endl
           << "Only a " << (xCurrent-x1start)*100/(x2end-x1start)
           << " % fraction of it was done.";
   G4Exception("G4OldMagIntDriver::WarnTooManyStep()", "GeomField1001",
               JustWarning, message);
}

// ---------------------------------------------------------

void
G4OldMagIntDriver::WarnEndPointTooFar (G4double endPointDist, 
                                     G4double h , 
                                     G4double eps,
                                     G4int    dbg)
{
  static G4ThreadLocal G4double maxRelError = 0.0;
  G4bool isNewMax, prNewMax;

  isNewMax = endPointDist > (1.0 + maxRelError) * h;
  prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
  if( isNewMax ) { maxRelError= endPointDist / h - 1.0; }

  if( dbg && (h > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) 
          && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+eps) ) ) )
  { 
    static G4ThreadLocal G4int noWarnings = 0;
    std::ostringstream message;
    if( (noWarnings++ < 10) || (dbg>2) )
    {
      message << "The integration produced an end-point which " << G4endl
              << "is further from the start-point than the curve length."
              << G4endl;
    }
    message << "  Distance of endpoints = " << endPointDist
            << ", curve length = " << h << G4endl
            << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
            << ", relative = " << (h-endPointDist) / h 
            << ", epsilon =  " << eps;
    G4Exception("G4OldMagIntDriver::WarnEndPointTooFar()", "GeomField1001",
                JustWarning, message);
  }
}

// ---------------------------------------------------------

void
G4OldMagIntDriver::OneGoodStep(      G4double y[],        // InOut
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
  G4double errmax_sq;
  G4double h, htemp, xnew ;

  G4double yerr[G4FieldTrack::ncompSVEC], ytemp[G4FieldTrack::ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  G4double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  G4double errpos_sq = 0.0;    // square of displacement error
  G4double errvel_sq = 0.0;    // square of momentum vector difference
  G4double errspin_sq = 0.0;   // square of spin vector difference

  const G4int max_trials=100; 

  G4ThreeVector Spin(y[9],y[10],y[11]);
  G4double spin_mag2 = Spin.mag2();
  G4bool hasSpin = (spin_mag2 > 0.0); 

  for (G4int iter=0; iter<max_trials; ++iter)
  {
    pIntStepper-> Stepper(y,dydx,h,ytemp,yerr); 
    //            *******
    G4double eps_pos = eps_rel_max * std::max(h, fMinimumStep); 
    G4double inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos); 

    // Evaluate accuracy
    //
    errpos_sq =  sqr(yerr[0]) + sqr(yerr[1]) + sqr(yerr[2]) ;
    errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

    // Accuracy for momentum
    G4double magvel_sq=  sqr(y[3]) + sqr(y[4]) + sqr(y[5]) ;
    G4double sumerr_sq =  sqr(yerr[3]) + sqr(yerr[4]) + sqr(yerr[5]) ; 
    if( magvel_sq > 0.0 )
    { 
       errvel_sq = sumerr_sq / magvel_sq; 
    }
    else
    {
       std::ostringstream message;
       message << "Found case of zero momentum." << G4endl
               << "- iteration= " << iter << "; h= " << h;
       G4Exception("G4OldMagIntDriver::OneGoodStep()",
                   "GeomField1001", JustWarning, message);
       errvel_sq = sumerr_sq; 
    }
    errvel_sq *= inv_eps_vel_sq;
    errmax_sq = std::max( errpos_sq, errvel_sq ); // Square of maximum error

    if( hasSpin )
    { 
      // Accuracy for spin
      errspin_sq =  ( sqr(yerr[9]) + sqr(yerr[10]) + sqr(yerr[11]) )
                    /  spin_mag2; // ( sqr(y[9]) + sqr(y[10]) + sqr(y[11]) );
      errspin_sq *= inv_eps_vel_sq;
      errmax_sq = std::max( errmax_sq, errspin_sq ); 
    }

    if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

    // Step failed; compute the size of retrial Step.
    htemp = GetSafety() * h * std::pow( errmax_sq, 0.5*GetPshrnk() );

    if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
    else  { h = 0.1*h; }                 // reduce stepsize, but no more
                                         // than a factor of 10
    xnew = x + h;
    if(xnew == x)
    {
      std::ostringstream message;
      message << "Stepsize underflow in Stepper !" << G4endl
              << "- Step's start x=" << x << " and end x= " << xnew 
              << " are equal !! " << G4endl
              << "  Due to step-size= " << h 
              << ". Note that input step was " << htry;
      G4Exception("G4OldMagIntDriver::OneGoodStep()",
                  "GeomField1001", JustWarning, message);
      break;
    }
  }

  // Compute size of next Step
  if (errmax_sq > errcon*errcon)
  { 
    hnext = GetSafety()*h*std::pow(errmax_sq, 0.5*GetPgrow());
  }
  else
  {
    hnext = max_stepping_increase*h ; // No more than a factor of 5 increase
  }
  x += (hdid = h);

  for(G4int k=0; k<fNoIntegrationVariables; ++k) { y[k] = ytemp[k]; }

  return;
}

//----------------------------------------------------------------------

// QuickAdvance just tries one Step - it does not ensure accuracy
//
G4bool G4OldMagIntDriver::QuickAdvance(G4FieldTrack& y_posvel,    // INOUT
                               const G4double      dydx[],  
                                     G4double      hstep,       // In
                                     G4double&     dchord_step,
                                     G4double&     dyerr_pos_sq,
                                     G4double&     dyerr_mom_rel_sq )  
{
  G4Exception("G4OldMagIntDriver::QuickAdvance()", "GeomField0001",
              FatalException, "Not yet implemented."); 

  // Use the parameters of this method, to please compiler
  //
  dchord_step = dyerr_pos_sq = hstep * hstep * dydx[0]; 
  dyerr_mom_rel_sq = y_posvel.GetPosition().mag2();
  return true;
}

//----------------------------------------------------------------------

G4bool G4OldMagIntDriver::QuickAdvance(G4FieldTrack& y_posvel,    // INOUT
                               const G4double      dydx[],  
                                     G4double      hstep,       // In
                                     G4double&     dchord_step,
                                     G4double&     dyerr )
{
  G4double dyerr_pos_sq, dyerr_mom_rel_sq;  
  G4double yerr_vec[G4FieldTrack::ncompSVEC],
           yarrin[G4FieldTrack::ncompSVEC], yarrout[G4FieldTrack::ncompSVEC]; 
  G4double s_start;
  G4double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

#ifdef  G4DEBUG_FIELD  
  G4FieldTrack startTrack( y_posvel );  // For debugging
#endif

  // Move data into array
  y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
  s_start = y_posvel.GetCurveLength();

  // Do an Integration Step
  pIntStepper-> Stepper(yarrin, dydx, hstep, yarrout, yerr_vec) ; 

  // Estimate curve-chord distance
  dchord_step= pIntStepper-> DistChord();

  // Put back the values.  yarrout ==> y_posvel
  y_posvel.LoadFromArray( yarrout, fNoIntegrationVariables );
  y_posvel.SetCurveLength( s_start + hstep );

#ifdef  G4DEBUG_FIELD
  if(fVerboseLevel>2)
  {
    G4cout << "G4MagIntDrv: Quick Advance" << G4endl;
    PrintStatus( yarrin, s_start, yarrout, s_start+hstep, hstep,  1); 
  }
#endif

  // A single measure of the error   
  //      TO-DO :  account for  energy,  spin, ... ? 
  vel_mag_sq   = ( sqr(yarrout[3])+sqr(yarrout[4])+sqr(yarrout[5]) );
  inv_vel_mag_sq = 1.0 / vel_mag_sq; 
  dyerr_pos_sq = ( sqr(yerr_vec[0])+sqr(yerr_vec[1])+sqr(yerr_vec[2]));
  dyerr_mom_sq = ( sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));
  dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

  // Calculate also the change in the momentum squared also ???
  // G4double veloc_square = y_posvel.GetVelocity().mag2();
  // ...

#ifdef RETURN_A_NEW_STEP_LENGTH
  // The following step cannot be done here because "eps" is not known.
  dyerr_len = std::sqrt( dyerr_len_sq ); 
  dyerr_len_sq /= eps ;

  // Look at the velocity deviation ?
  //  sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));

  // Set suggested new step
  hstep = ComputeNewStepSize( dyerr_len, hstep);
#endif

  if( dyerr_pos_sq > ( dyerr_mom_rel_sq * sqr(hstep) ) )
  {
    dyerr = std::sqrt(dyerr_pos_sq);
  }
  else
  {
    // Scale it to the current step size - for now
    dyerr = std::sqrt(dyerr_mom_rel_sq) * hstep;
  }

#ifdef  G4DEBUG_FIELD  
  // For debugging
  G4cout // << "G4MagInt_Driver::"
         << "QuickAdvance" << G4endl
         << " Input:  hstep= " << hstep       << G4endl
         << "         track= " << startTrack  << G4endl
         << " Output: track= " << y_posvel    << G4endl
         << "         d_chord = " << dchord_step
         << " dyerr = " << dyerr << G4endl;
#endif

  return true;
}

// --------------------------------------------------------------------------

#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
G4bool  G4OldMagIntDriver::QuickAdvance(G4double  yarrin[],    // In
                                const G4double  dydx[],  
                                      G4double  hstep,       // In
                                      G4double  yarrout[],
                                      G4double& dchord_step,
                                      G4double& dyerr )      // In length
{
  G4Exception("G4OldMagIntDriver::QuickAdvance()", "GeomField0001",
              FatalException, "Not yet implemented.");
  dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
  yarrout[0]= yarrin[0];
}
#endif 

// --------------------------------------------------------------------------

// This method computes new step sizes - but does not limit changes to
// within  certain factors
// 
G4double G4OldMagIntDriver::
ComputeNewStepSize(G4double  errMaxNorm,    // max error  (normalised)
                   G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if(errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPshrnk()) ;
  }
  else if(errMaxNorm > 0.0 )
  {
    // Compute size of next Step for a successful step
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPgrow()) ;
  }
  else
  {
    // if error estimate is zero (possible) or negative (dubious)
    hnew = max_stepping_increase * hstepCurrent; 
  }

  return hnew;
}

// ---------------------------------------------------------------------------

// This method computes new step sizes limiting changes within certain factors
// 
// It shares its logic with AccurateAdvance.
// They are kept separate currently for optimisation.
//
G4double 
G4OldMagIntDriver::ComputeNewStepSize_WithinLimits( 
                          G4double  errMaxNorm,    // max error  (normalised)
                          G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if (errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPshrnk()) ;
  
    if (hnew < max_stepping_decrease*hstepCurrent)
    {
      hnew = max_stepping_decrease*hstepCurrent ;
                         // reduce stepsize, but no more
                         // than this factor (value= 1/10)
    }
  }
  else
  {
    // Compute size of next Step for a successful step
    if (errMaxNorm > errcon)
     { hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPgrow()); }
    else  // No more than a factor of 5 increase
     { hnew = max_stepping_increase * hstepCurrent; }
  }
  return hnew;
}

// ---------------------------------------------------------------------------

void G4OldMagIntDriver::PrintStatus( const G4double* StartArr,  
                                         G4double  xstart,
                                   const G4double* CurrentArr, 
                                         G4double  xcurrent,
                                         G4double  requestStep, 
                                         G4int     subStepNo )
  // Potentially add as arguments:  
  //                                 <dydx>           - as Initial Force
  //                                 stepTaken(hdid)  - last step taken
  //                                 nextStep (hnext) - proposal for size
{
   G4FieldTrack  StartFT(G4ThreeVector(0,0,0),
                 G4ThreeVector(0,0,0), 0., 0., 0., 0. );
   G4FieldTrack  CurrentFT (StartFT);

   StartFT.LoadFromArray( StartArr, fNoIntegrationVariables); 
   StartFT.SetCurveLength( xstart);
   CurrentFT.LoadFromArray( CurrentArr, fNoIntegrationVariables); 
   CurrentFT.SetCurveLength( xcurrent );

   PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
}

// ---------------------------------------------------------------------------

void G4OldMagIntDriver::PrintStatus(const G4FieldTrack& StartFT,
                                  const G4FieldTrack& CurrentFT, 
                                        G4double      requestStep, 
                                        G4int         subStepNo)
{
    G4int verboseLevel= fVerboseLevel;
    const G4int noPrecision = 5;
    G4long oldPrec= G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);

    const G4ThreeVector StartPosition=       StartFT.GetPosition();
    const G4ThreeVector StartUnitVelocity=   StartFT.GetMomentumDir();
    const G4ThreeVector CurrentPosition=     CurrentFT.GetPosition();
    const G4ThreeVector CurrentUnitVelocity= CurrentFT.GetMomentumDir();

    G4double  DotStartCurrentVeloc= StartUnitVelocity.dot(CurrentUnitVelocity);

    G4double step_len= CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
    G4double subStepSize = step_len;
     
    if( (subStepNo <= 1) || (verboseLevel > 3) )
    {
       subStepNo = - subStepNo;        // To allow printing banner

       G4cout << std::setw( 6)  << " " << std::setw( 25)
              << " G4OldMagIntDriver: Current Position  and  Direction" << " "
              << G4endl; 
       G4cout << std::setw( 5) << "Step#" << " "
              << std::setw( 7) << "s-curve" << " "
              << std::setw( 9) << "X(mm)" << " "
              << std::setw( 9) << "Y(mm)" << " "  
              << std::setw( 9) << "Z(mm)" << " "
              << std::setw( 8) << " N_x " << " "
              << std::setw( 8) << " N_y " << " "
              << std::setw( 8) << " N_z " << " "
              << std::setw( 8) << " N^2-1 " << " "
              << std::setw(10) << " N(0).N " << " "
              << std::setw( 7) << "KinEner " << " "
              << std::setw(12) << "Track-l" << " "   // Add the Sub-step ??
              << std::setw(12) << "Step-len" << " " 
              << std::setw(12) << "Step-len" << " " 
              << std::setw( 9) << "ReqStep" << " "  
              << G4endl;
    }

    if( (subStepNo <= 0) )
    {
      PrintStat_Aux( StartFT,  requestStep, 0., 
                       0,        0.0,         1.0);
    }

    if( verboseLevel <= 3 )
    {
      G4cout.precision(noPrecision);
      PrintStat_Aux( CurrentFT, requestStep, step_len, 
                     subStepNo, subStepSize, DotStartCurrentVeloc );
    }

    G4cout.precision(oldPrec);
}

// ---------------------------------------------------------------------------

void G4OldMagIntDriver::PrintStat_Aux(const G4FieldTrack& aFieldTrack,
                                          G4double      requestStep, 
                                          G4double      step_len,
                                          G4int         subStepNo,
                                          G4double      subStepSize,
                                          G4double      dotVeloc_StartCurr)
{
    const G4ThreeVector Position = aFieldTrack.GetPosition();
    const G4ThreeVector UnitVelocity = aFieldTrack.GetMomentumDir();
 
    if( subStepNo >= 0)
    {
       G4cout << std::setw( 5) << subStepNo << " ";
    }
    else
    {
       G4cout << std::setw( 5) << "Start" << " ";
    }
    G4double curveLen= aFieldTrack.GetCurveLength();
    G4cout << std::setw( 7) << curveLen;
    G4cout << std::setw( 9) << Position.x() << " "
           << std::setw( 9) << Position.y() << " "
           << std::setw( 9) << Position.z() << " "
           << std::setw( 8) << UnitVelocity.x() << " "
           << std::setw( 8) << UnitVelocity.y() << " "
           << std::setw( 8) << UnitVelocity.z() << " ";
    G4long oldprec= G4cout.precision(3);
    G4cout << std::setw( 8) << UnitVelocity.mag2()-1.0 << " ";
    G4cout.precision(6);
    G4cout << std::setw(10) << dotVeloc_StartCurr << " ";
    G4cout.precision(oldprec);
    G4cout << std::setw( 7) << aFieldTrack.GetKineticEnergy();
    G4cout << std::setw(12) << step_len << " ";

    static G4ThreadLocal G4double oldCurveLength = 0.0;
    static G4ThreadLocal G4double oldSubStepLength = 0.0;
    static G4ThreadLocal G4int oldSubStepNo = -1;

    G4double subStep_len = 0.0;
    if( curveLen > oldCurveLength )
    {
      subStep_len= curveLen - oldCurveLength;
    }
    else if (subStepNo == oldSubStepNo)
    {
      subStep_len= oldSubStepLength;
    }
    oldCurveLength= curveLen;
    oldSubStepLength= subStep_len;

    G4cout << std::setw(12) << subStep_len << " "; 
    G4cout << std::setw(12) << subStepSize << " "; 
    if( requestStep != -1.0 )
    {
      G4cout << std::setw( 9) << requestStep << " ";
    }
    else
    {
       G4cout << std::setw( 9) << " InitialStep " << " ";
    }
    G4cout << G4endl;
}

// ---------------------------------------------------------------------------

void G4OldMagIntDriver::PrintStatisticsReport()
{
  G4int noPrecBig = 6;
  G4long oldPrec = G4cout.precision(noPrecBig);

  G4cout << "G4OldMagIntDriver Statistics of steps undertaken. " << G4endl;
  G4cout << "G4OldMagIntDriver: Number of Steps: "
         << " Total= " <<  fNoTotalSteps
         << " Bad= "   <<  fNoBadSteps 
         << " Small= " <<  fNoSmallSteps 
         << " Non-initial small= " << (fNoSmallSteps-fNoInitialSmallSteps)
         << G4endl;
 G4cout.precision(oldPrec);
}
 
// ---------------------------------------------------------------------------

void G4OldMagIntDriver::SetSmallestFraction(G4double newFraction)
{
  if( (newFraction > 1.e-16) && (newFraction < 1e-8) )
  {
    fSmallestFraction= newFraction;
  }
  else
  { 
    std::ostringstream message;
    message << "Smallest Fraction not changed. " << G4endl
            << "  Proposed value was " << newFraction << G4endl
            << "  Value must be between 1.e-8 and 1.e-16";
    G4Exception("G4OldMagIntDriver::SetSmallestFraction()",
                "GeomField1001", JustWarning, message);
  }
}

void G4OldMagIntDriver::
GetDerivatives(const G4FieldTrack& y_curr, G4double* dydx) const
{
    G4double ytemp[G4FieldTrack::ncompSVEC];
    y_curr.DumpToArray(ytemp);
    pIntStepper->RightHandSide(ytemp, dydx);
      // Avoid virtual call for GetStepper
      // Was: GetStepper()->ComputeRightHandSide(ytemp, dydx);
}

void G4OldMagIntDriver::GetDerivatives(const G4FieldTrack& track,
                                     G4double dydx[],
                                     G4double field[]) const
{
    G4double ytemp[G4FieldTrack::ncompSVEC];
    track.DumpToArray(ytemp);
    pIntStepper->RightHandSide(ytemp, dydx, field);
}

G4EquationOfMotion* G4OldMagIntDriver::GetEquationOfMotion()
{
    return pIntStepper->GetEquationOfMotion();
}

void G4OldMagIntDriver::SetEquationOfMotion(G4EquationOfMotion *equation)
{
    pIntStepper->SetEquationOfMotion(equation);
}

const G4MagIntegratorStepper* G4OldMagIntDriver::GetStepper() const
{
    return pIntStepper;
}

G4MagIntegratorStepper* G4OldMagIntDriver::GetStepper()
{
    return pIntStepper;
}

void G4OldMagIntDriver::
RenewStepperAndAdjust(G4MagIntegratorStepper* pItsStepper)
{  
    pIntStepper = pItsStepper; 
    ReSetParameters();
}

void G4OldMagIntDriver::StreamInfo( std::ostream& os ) const
{
    os << "State of G4OldMagIntDriver: " << std::endl;
    os << "  Max number of Steps = " << fMaxNoSteps 
       << "    (base # = " << fMaxStepBase << " )" << std::endl;
    os << "  Safety factor       = " << safety  << std::endl;
    os << "  Power - shrink      = " << pshrnk << std::endl;
    os << "  Power - grow        = " << pgrow << std::endl;
    os << "  threshold (errcon)  = " << errcon << std::endl;

    os << "    fMinimumStep =      " << fMinimumStep  << std::endl;
    os << "    Smallest Fraction = " << fSmallestFraction << std::endl;

    os << "    No Integrat Vars  = " << fNoIntegrationVariables << std::endl;
    os << "    Min No Vars       = " << fMinNoVars << std::endl;
    os << "    Num-Vars          = " << fNoVars << std::endl;
    
    os << "    verbose level     = " << fVerboseLevel << std::endl;

    // auto p= const_cast<G4OldMagIntDriver*>(this);
    G4bool does = // p->DoesReIntegrate();
          const_cast<G4OldMagIntDriver*>(this)->DoesReIntegrate();    
    os << "    Reintegrates      = " << does  << std::endl;
}
