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
// $Id: G4MagErrorStepper.cc,v 1.9 2001-07-11 09:59:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4MagErrorStepper.hh"
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

G4MagErrorStepper::~G4MagErrorStepper()
{
   delete[] yMiddle;
   delete[] dydxMid;
   delete[] yInitial;
   delete[] yOneStep;
}

void
G4MagErrorStepper::Stepper( const G4double yInput[],
		            const G4double dydx[],
		                  G4double hstep,
		                  G4double yOutput[],
		                  G4double yError []      )
{  
   const G4int nvar = this->GetNumberOfVariables() ;

   G4int i;
   // correction for Richardson Extrapolation.
   G4double  correction = 1. / ( (1 << IntegratorOrder()) -1 );
   
   //  Saving yInput because yInput and yOutput can be aliases for same array

   for(i=0;i<nvar;i++) yInitial[i]=yInput[i];

   G4double halfStep = hstep * 0.5; 

   // Do two half steps

   DumbStepper  (yInitial,  dydx,   halfStep, yMiddle);
   RightHandSide(yMiddle, dydxMid);    
   DumbStepper  (yMiddle, dydxMid, halfStep, yOutput); 

   // Store midpoint, chord calculation

   fMidPoint = G4ThreeVector( yMiddle[0],  yMiddle[1],  yMiddle[2]); 

   // Do a full Step
   DumbStepper(yInitial, dydx, hstep, yOneStep);
   for(i=0;i<nvar;i++) {
      yError [i] = yOutput[i] - yOneStep[i] ;
      yOutput[i] += yError[i]*correction ;  // Provides accuracy increased
                                            // by 1 order via the 
                                            // Richardson Extrapolation  
   }

   fInitialPoint = G4ThreeVector( yInitial[0], yInitial[1], yInitial[2]); 
   fFinalPoint   = G4ThreeVector( yOutput[0],  yOutput[1],  yOutput[2]); 

   return ;
}



G4double
G4MagErrorStepper::DistChord() const 
{
  // Estimate the maximum distance from the curve to the chord
  //
  //  We estimate this using the distance of the midpoint to 
  //  chord (the line between 
  // 
  //  Method below is good only for angle deviations < 2 pi, 
  //   This restriction should not a problem for the Runge cutta methods, 
  //   which generally cannot integrate accurately for large angle deviations.
  G4double distLine, distChord; 

  if (fInitialPoint != fFinalPoint) {
     distLine= G4LineSection::Distline( fMidPoint, fInitialPoint, fFinalPoint );
     // This is a class method that gives distance of Mid 
     //  from the Chord between the Initial and Final points.

     distChord = distLine;
  }else{
     distChord = (fMidPoint-fInitialPoint).mag();
  }

  return distChord;
}
 
