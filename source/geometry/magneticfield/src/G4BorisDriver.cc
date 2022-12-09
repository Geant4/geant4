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
// G4BorisDriver
//
// Class description:
//
//   G4BorisDriver is a driver class using the second order Boris 
// method to integrate the equation of motion.
// 
//
// Author: Divyansh Tiwari, Google Summer of Code 2022
// Supervision: John Apostolakis,Renee Fatemi, Soon Yung Jun 
// --------------------------------------------------------------------
#include <cassert>

#include "G4BorisDriver.hh"

#include "G4SystemOfUnits.hh"
#include "G4LineSection.hh"
#include "G4FieldUtils.hh"

const G4double G4BorisDriver::fErrorConstraintShrink = std::pow(
            fMaxSteppingDecrease / fSafetyFactor, 1. / fPowerShrink);

const G4double G4BorisDriver::fErrorConstraintGrow = std::pow(
            fMaxSteppingIncrease / fSafetyFactor, 1. / fPowerGrow);

// --------------------------------------------------------------------------

G4BorisDriver::
G4BorisDriver( G4double hminimum, G4BorisScheme* Boris,
               G4int numberOfComponents, bool verbosity )
  : fMinimumStep(hminimum),
    fVerbosity(verbosity),
    boris(Boris)
    // , interval_sequence{2,4}
{    
    assert(boris->GetNumberOfVariables() == numberOfComponents);
    
    if(boris->GetNumberOfVariables() != numberOfComponents)
    {
      std::ostringstream msg;
      msg << "Disagreement in number of variables = "
          << boris->GetNumberOfVariables()
          << " vs no of components = " << numberOfComponents;
      G4Exception("G4BorisDriver Constructor:",
                  "GeomField1001", FatalException, msg);       
    }

}

// --------------------------------------------------------------------------

G4bool G4BorisDriver::AccurateAdvance( G4FieldTrack& track,
                                       G4double hstep,
                                       G4double  epsilon,
                                       G4double  hinitial )
{
   // Specification: Driver with adaptive stepsize control.
   // Integrate starting values at y_current over hstep x2 with (relative) accuracy 'eps'.
   // On output 'track' is replaced by values at the end of the integration interval. 
   
   //  Ensure that hstep > 0
   if(hstep == 0)
   {
      std::ostringstream message;
      message << "Proposed step is zero; hstep = " << hstep << " !";
      G4Exception("G4BorisDriver::AccurateAdvance()",
                  "GeomField1001", JustWarning, message);
      return true;
   }
   if(hstep < 0)
   {
      std::ostringstream message;
      message << "Invalid run condition." << G4endl
              << "Proposed step is negative; hstep = " << hstep << G4endl
              << "Requested step cannot be negative! Aborting event.";
      G4Exception("G4BorisDriver::AccurateAdvance()",
                  "GeomField0003", EventMustBeAborted, message);
      return false;
   }

   if( hinitial == 0.0 ) { hinitial = hstep; }   
   if( hinitial < 0.0 ) { hinitial = std::fabs( hinitial ); }
   // G4double htrial = std::min( hstep, hinitial );
   G4double htrial = hstep;
   // Decide first step size      

   // G4int noOfSteps = h/hstep;
   
    // integration variables
    //
   track.DumpToArray(yCurrent);

   const G4double restMass = track.GetRestMass();
   const G4double charge = track.GetCharge()*e_SI;
   const G4int    nvar= GetNumberOfVariables();
   
   // copy non-integration variables to out array
   //
   std::memcpy(yOut + nvar,
               yCurrent + nvar,
               sizeof(G4double)*(G4FieldTrack::ncompSVEC-nvar));
   
   G4double curveLength = track.GetCurveLength();  // starting value
   const G4double endCurveLength = curveLength + hstep;

   // -- Initial version: Did it in one step -- did not account for errors !!!
   // G4FieldTrack yFldTrk(track);
   // yFldTrk.LoadFromArray(yCurrent, G4FieldTrack::ncompSVEC);
   // yFldTrk.SetCurveLength(curveLength);
   // G4double dchord_step, dyerr_len;   
   // QuickAdvance(yFldTrk, dydxCurrent, htrial, dchord_step, dyerr_len);

   const G4double hThreshold = 
        std::max(epsilon * hstep, fSmallestFraction * curveLength);

   G4double htry= htrial;
   
   for (G4int nstp = 0; nstp < fMaxNoSteps; ++nstp)
   {
      G4double hdid= 0.0, hnext=0.0;
        
      OneGoodStep(yCurrent, curveLength, htry, epsilon, restMass, charge, hdid, hnext);
      //*********

      // Simple check: move (distance of displacement) is smaller than length along curve!
      const G4ThreeVector StartPos = field_utils::makeVector(yCurrent, field_utils::Value3D::Position);      
      const G4ThreeVector EndPos = field_utils::makeVector(yCurrent, field_utils::Value3D::Position);
      CheckStep(EndPos, StartPos, hdid);
      
      // Check   1. for finish   and    2. *avoid* numerous small last steps
      if (curveLength >= endCurveLength || htry < hThreshold)
      {
         break; 
      }
      
      htry = std::max(hnext, fMinimumStep);
      if (curveLength + htry > endCurveLength)
      {
         htry = endCurveLength - curveLength;
      }
   }

   // upload new state   
   track.LoadFromArray(yCurrent, G4FieldTrack::ncompSVEC);
   track.SetCurveLength(curveLength);

   return true;
}

// --------------------------------------------------------------------------

void G4BorisDriver::OneGoodStep(G4double y[],           // InOut
                                G4double& curveLength,  // InOut
                                G4double  htry,
                                G4double  epsilon_rel,
                                G4double  restMass,
                                G4double  charge,
                                G4double& hdid,        // Out
                                G4double& hnext)       // Out
{
    G4double error2 = DBL_MAX;
    G4double yerr[G4FieldTrack::ncompSVEC], ytemp[G4FieldTrack::ncompSVEC];

    G4double h = htry;

    const G4int max_trials = 100; 

    for (G4int iter = 0; iter < max_trials; ++iter)
    {
        boris->StepWithErrorEstimate(y, restMass, charge, h, ytemp, yerr);
        
        error2 = field_utils::relativeError2(y, yerr, std::max(h, fMinimumStep),
                                             epsilon_rel);
        if (error2 <= 1.0)
        {
            break; 
        }

        h = ShrinkStepSize2(h, error2);

        G4double xnew = curveLength + h;
        if(xnew == curveLength)
        {
            std::ostringstream message;
            message << "Stepsize underflow in Stepper !" << G4endl
                    << "- Step's start x=" << curveLength
                    << " and end x= " << xnew 
                    << " are equal !! " << G4endl
                    << "  Due to step-size= " << h 
                    << ". Note that input step was " << htry;
            G4Exception("G4IntegrationDriver::OneGoodStep()",
                        "GeomField1001", JustWarning, message);
            break;
        }
    }

    hnext = GrowStepSize2(h, error2);
    curveLength += (hdid = h);

    field_utils::copy(y, ytemp, GetNumberOfVariables());
}

// ===========------------------------------------------------------===========

G4bool G4BorisDriver::
QuickAdvance( G4FieldTrack& track, const G4double /*dydx*/[],
              G4double hstep, G4double& missDist, G4double& dyerr)
{
    const auto nvar = boris->GetNumberOfVariables();

    track.DumpToArray(yIn);
    const G4double curveLength = track.GetCurveLength();

    // call the boris method for step length hstep
    G4double restMass = track.GetRestMass();
    G4double charge = track.GetCharge()*e_SI;

    // boris->DoStep(restMass, charge, yIn,  yMid, hstep*0.5);
    // boris->DoStep(restMass, charge, yMid, yOut, hstep*0.5);  // Use mid-point !!
    boris->StepWithMidAndErrorEstimate(yIn, restMass, charge, hstep,
                                       yMid, yOut, yError);
      // Same, and also return mid-point evaluation
    
    // How to calculate chord length??
    const auto mid = field_utils::makeVector(yMid,
                     field_utils::Value3D::Position);
    const auto in  = field_utils::makeVector(yIn,
                     field_utils::Value3D::Position);
    const auto out = field_utils::makeVector(yOut,
                     field_utils::Value3D::Position);

    missDist = G4LineSection::Distline(mid, in, out);

    dyerr = field_utils::absoluteError(yOut, yError, hstep);

    // copy non-integrated variables to output array
    //
    std::memcpy(yOut + nvar, yIn + nvar,
                sizeof(G4double) * (G4FieldTrack::ncompSVEC - nvar));

    // set new state
    //
    track.LoadFromArray(yOut, G4FieldTrack::ncompSVEC);
    track.SetCurveLength(curveLength +  hstep);

    return true;
}

// --------------------------------------------------------------------------------

void G4BorisDriver::
GetDerivatives( const G4FieldTrack& yTrack, G4double dydx[]) const
{
   G4double  EBfieldValue[6];
   GetDerivatives(yTrack, dydx, EBfieldValue);
}

// --------------------------------------------------------------------------------

void G4BorisDriver::
GetDerivatives( const G4FieldTrack& yTrack, G4double dydx[],
                G4double EBfieldValue[]) const
{
   // G4Exception("G4BorisDriver::GetDerivatives()",
   //             "GeomField0003", FatalException, "This method is not implemented.");
   G4double ytemp[G4FieldTrack::ncompSVEC];
   yTrack.DumpToArray(ytemp);
   GetEquationOfMotion()->EvaluateRhsReturnB(ytemp, dydx, EBfieldValue);
}

// --------------------------------------------------------------------------------

G4double G4BorisDriver::ShrinkStepSize2(G4double h, G4double error2) const
{
   if (error2 > fErrorConstraintShrink * fErrorConstraintShrink)
   {
      return fMaxSteppingDecrease * h;
   }
   return fSafetyFactor * h * std::pow(error2, 0.5 * fPowerShrink);
}

// --------------------------------------------------------------------------------

G4double G4BorisDriver::GrowStepSize2(G4double h, G4double error2) const
// Given the square of the relative error, 
{
    if (error2 < fErrorConstraintGrow * fErrorConstraintGrow)
    {
        return fMaxSteppingIncrease * h;
    }
    return fSafetyFactor * h * std::pow(error2, 0.5 * fPowerGrow);
}

// --------------------------------------------------------------------------------

void G4BorisDriver::SetEquationOfMotion(G4EquationOfMotion* /*equation*/ )
{
   G4Exception("G4BorisDriver::SetEquationOfMotion()", "GeomField0003", FatalException,
               "This method is not implemented. BorisDriver/Stepper should keep its equation");
}

// --------------------------------------------------------------------------------

void
G4BorisDriver::StreamInfo( std::ostream& os ) const
{
  os << "State of G4BorisDriver: " << std::endl;
  os << "   Method is implemented, but gives no information. " << std::endl;
}
