// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagErrorStepper.cc,v 1.5 1999-04-19 17:20:30 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4MagErrorStepper.hh"
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

void
G4MagErrorStepper::Stepper( const G4double yInput[],
		     const G4double dydx[],
		     const G4double hstep,
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
G4MagErrorStepper::DistChord()   const 
{
  // Soon: must check whether h/R > 2 pi  !!
  //  Method below is good only for < 2 pi

  return G4LineSection::Distline( fMidPoint, fInitialPoint, fFinalPoint );
  // This is a class method that gives distance of Mid 
  //  from the Chord between the Initial and Final points.
}
 
