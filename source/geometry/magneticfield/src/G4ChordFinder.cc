// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ChordFinder.cc,v 1.7 1999-07-19 17:31:44 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// 25.02.97 John Apostolakis,  design and implimentation 
// 05.03.97 V. Grichine , style modification

#include "G4ChordFinder.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
// #include "G4Field.hh"
                                       // #include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"

//  For the moment fDeltaChord is a constant!

const G4double G4ChordFinder::fDefaultDeltaChord  = 3. * mm; 

// ..........................................................................

G4ChordFinder::G4ChordFinder( G4MagneticField*        theMagField,
		              G4double                stepMinimum, 
		              G4MagIntegratorStepper* pItsStepper ) // A default one
     : fDeltaChord( fDefaultDeltaChord )
{
  //  Construct the Chord Finder
  //  by creating in inverse order the  Driver, the Stepper and EqRhs ...
  // G4Mag_EqRhs *
  fEquation = new G4Mag_UsualEqRhs(theMagField); // Should move q, p to 
                                                     //G4FieldTrack ??
  // --->>  Charge    Q = 0 
  // --->>  Momentum  P = 1       NOMINAL VALUES !!!!!!!!!!!!!!!!!!

  if( pItsStepper == 0 )
  { 
     pItsStepper = fDriversStepper = new G4ClassicalRK4(fEquation);
     fAllocatedStepper= true;
  }
  else
  {
     fAllocatedStepper= false; 
  }
  fIntgrDriver = new G4MagInt_Driver(stepMinimum, 
				     pItsStepper, 
				     pItsStepper->GetNumberOfVariables() );
}

// ......................................................................

G4ChordFinder::~G4ChordFinder()
{
  delete   fEquation; // fIntgrDriver->pIntStepper->theEquation_Rhs;
  if( fAllocatedStepper)
  { 
     delete fDriversStepper; 
  }                                //  fIntgrDriver->pIntStepper;}
  delete   fIntgrDriver; 
}

// ......................................................................

G4double 
G4ChordFinder::AdvanceChordLimited(   G4FieldTrack& yCurrent,
				   const  G4double     stepMax,
				   const  G4double     epsStep )
{
  G4double stepPossible;
  G4double dyErr;
  G4FieldTrack yEnd( yCurrent);

  G4bool dbg= false; 

#ifdef G4VERBOSE
  if( dbg ) 
    G4cerr << "Entered FindNextChord Limited with:\n yCurrent: " << yCurrent
	   << " and initial Step=stepMax=" <<  stepMax << " mm. " << endl;
#endif

  stepPossible= FindNextChord(yCurrent, stepMax, yEnd, dyErr, epsStep);
  G4bool good_advance;
  if ( dyErr < epsStep * stepPossible )
  {
     // Accept this accuracy.
     yCurrent = yEnd;
     good_advance = true; 
  }
  else
  {
     // Advance more accurately to "end of chord"
     good_advance = fIntgrDriver->AccurateAdvance(yCurrent, stepPossible, epsStep);
     #ifdef G4VERBOSE
     if (dbg) G4cerr << "Accurate advance to end of chord attemped"
		       << "with result " << good_advance << endl ;
     #endif
     if ( ! good_advance ){ 
       // In this case the driver could not do the full distance
       stepPossible= yCurrent.GetCurveLength();
     }
  }

#ifdef G4VERBOSE
  if( dbg ) G4cerr << "Exiting FindNextChord Limited with:\n yCurrent: " 
		 << yCurrent<< endl; 
#endif

  return stepPossible;
}

// ..............................................................................

G4double
G4ChordFinder::FindNextChord( const  G4FieldTrack  yStart,
	                      const  G4double     stepMax,
	                      G4FieldTrack&   yEnd,      //  Endpoint
	                      G4double&      dyErr,      //  Error of endpoint 
			      G4double     epsStep )
	    
// Returns Length of Step taken
{
  // G4int       stepRKnumber=0;
  G4FieldTrack yCurrent=  yStart;  
  G4double    stepTrial= stepMax;
  G4double    dydx[G4FieldTrack::ncompSVEC]; 

  //  1.)  Try to "leap" to end of interval
  //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
  //     2a.)  If d_chord is not good enough, find one that is.
  
  G4bool    validEndPoint= false,  dbg= false;
  G4double  dChordStep;

  fIntgrDriver-> GetDerivatives( yCurrent, dydx )  ;

  do
  { 
     yCurrent = yStart;    // Always start from initial point

     fIntgrDriver->QuickAdvance( yCurrent, dydx, stepTrial, dChordStep, dyErr);

#ifdef G4VERBOSE
     if( dbg ) {
        G4cerr << "Returned from QuickAdvance with: yCur=" << yCurrent << endl;
        G4cerr << " dChordStep= "<< dChordStep <<" dyErr=" << dyErr << endl; 
     }
#endif

     // We check whether the criterion is met here.
     validEndPoint = AcceptableMissDist(dChordStep); 
                      //  && (dyErr < eps) ;

     if( ! validEndPoint ) {
         // This is needed to decide new step size until QuickAdvance does it
	 stepTrial = NewStep(stepTrial, dChordStep );

	 // Get the driver to calculate the new step size, if it is needed
	 // stepTrial= fIntgrDriver->ComputeNewStepSize( dyErr/epsStep, stepTrial);
#ifdef G4VERBOSE
	 if( dbg ) 
	   G4cerr << "Dchord too big. Trying new hstep=" << stepTrial << endl;
#endif
     }
 
  }
  while( ! validEndPoint );   // End of do-while  RKD 

  yEnd=  yCurrent;  
  return stepTrial; 
}

// ...........................................................................

G4double G4ChordFinder::NewStep( 
		      const G4double stepTrialOld, 
		      const G4double dChordStep  )  // Current dchord achieved.
		   
{
  G4double stepTrial;

  if ( dChordStep > 1000. * fDeltaChord ){
        stepTrial= stepTrialOld * 0.03;   
  }else{
     if ( dChordStep > 100. * fDeltaChord ){
	stepTrial= stepTrialOld * 0.1;   
     }else{
        // Keep halving the length until dChordStep OK
	stepTrial= stepTrialOld * 0.5;   
     }
  }

  // A more sophisticated chord-finder could figure out a better
  //   stepTrial, from dChordStep and the required d_geometry
  //   eg
  //      Calculate R, r_helix (eg at orig point)
  //      if( stepTrial < 2 pi  R )
  //          stepTrial = R arc_cos( 1 - fDeltaChord / r_helix )
  //      else    
  //          ??

  return stepTrial;
}

//
//   Given a starting curve point A (CurveA_PointVelocity),  a later 
//  curve point B (CurveB_PointVelocity) and a point E which is (generally)
//  not on the curve, find and return a point F which is on the curve and 
//  which is close to E. While advancing towards F utilise eps_step 
//  as a measure of the relative accuracy of each Step.
  
G4FieldTrack G4ChordFinder::ApproxCurvePointV( 
			      const G4FieldTrack& CurveA_PointVelocity, 
			      const G4FieldTrack& CurveB_PointVelocity, 
			      const G4ThreeVector& CurrentE_Point,
			      const G4double eps_step)
{
  // 1st implementation:
  //    if r=|AE|/|AB|, and s=true path lenght (AB)
  //    return the point that is r*s along the curve!

  G4FieldTrack    Current_PointVelocity= CurveA_PointVelocity; 

  G4ThreeVector  CurveA_Point= CurveA_PointVelocity.Position();
  G4ThreeVector  CurveB_Point= CurveB_PointVelocity.Position();

  G4ThreeVector  ChordAB_Vector= CurveB_Point   - CurveA_Point;
  G4ThreeVector  ChordAE_Vector= CurrentE_Point - CurveA_Point;

  G4double       ABdist= ChordAB_Vector.mag();
  G4double  curve_length;  //  A curve length  of AB
  G4double  AE_fraction; 
  
  curve_length= 
       CurveB_PointVelocity.CurveS() - CurveA_PointVelocity.CurveS();  

  // const 
  G4double  integrationInaccuracyLimit= max( perMillion, 0.5*eps_step ); 
  if( curve_length < ABdist * (1. - integrationInaccuracyLimit) ){ 
//  #ifdef G4DEBUG
    G4cerr << " Warning in G4ChordFinder::ApproxCurvePoint: " << endl <<
      " The two points are further apart than the curve length " << endl <<
      " Dist = "         << ABdist  << 
      " curve length = " << curve_length 
	   << " relativeDiff = " << (curve_length-ABdist)/ABdist 
	   << endl;
//  #endif
    if( curve_length < ABdist * (1. - 10*eps_step) ) {
//    #ifdef G4DEBUG
      G4cerr << " ERROR: the size of the above difference exceeds allowed limits.  Aborting." 
	     << endl;
//    #endif
      G4Exception("G4ChordFinder::ApproxCurvePoint> Unphysical curve length.");
    }
    // Take default corrective action: 
    //    -->  adjust the maximum curve length. 
    //  NOTE: this case only happens for relatively straight paths.
    curve_length = ABdist; 
  }

  G4double  new_st_length; 

  if ( ABdist > 0.0 ){
     AE_fraction = ChordAE_Vector.mag() / ABdist;
  }else{
     G4cerr << " Error in G4ChordFinder::ApproxCurvePoint: A and B are the same point\n" <<
      " Chord AB length = " << ChordAE_Vector.mag()  << endl << endl;
     AE_fraction = 0.5;                         // Guess .. ?; 
  }
  
  if( (AE_fraction> 1.0 + perMillion) || (AE_fraction< 0.) ){
    G4cerr << " G4ChordFinder::ApproxCurvePointV: Warning: Anomalous condition:AE > AB or AE/AB <= 0 " << endl <<
      "   AE_fraction = " <<  AE_fraction << endl <<
      "   Chord AE length = " << ChordAE_Vector.mag()  << endl << 
      "   Chord AB length = " << ABdist << endl << endl;
    G4cerr << " OK if this condition occurs after a recalculation of 'B'" << endl
	   << " Otherwise it is an error. " << endl ; 
     // This course can now result if B has been re-evaluated, 
     //   without E being recomputed   (1 July 99)
     //  In this case this is not a "real error" - but it undesired
     //   and we cope with it by a default corrective action ...
     AE_fraction = 0.5;                         // Default value
  }

  new_st_length= AE_fraction * curve_length; 

  G4bool good_advance;
  if ( AE_fraction > 0.0 ) { 
     good_advance = 
      fIntgrDriver->AccurateAdvance(Current_PointVelocity, 
				    new_st_length,
				    eps_step ); // Relative accuracy
     // In this case it does not matter if it cannot advance the full distance
  }

  // If there was a memory of the step_length actually require at the start 
  // of the integration Step, this could be re-used ...

  return Current_PointVelocity;
}


