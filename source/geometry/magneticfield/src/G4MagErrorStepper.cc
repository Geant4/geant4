// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagErrorStepper.cc,v 1.1 1999-01-07 16:07:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4MagErrorStepper.hh"
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

void
G4MagErrorStepper::Stepper( const G4double yInput[],
		     const G4double dydx[],
		     const G4double hstep,
		     G4double yOut[],
		     G4double yErr[]      )
{  
   const G4int nvar = 6 ;

   G4int i;
   // correction for Richardson Extrapolation.
   G4double  correction = 1. / ( (1 << IntegratorOrder()) -1 );
   
   G4double yTemp[7], dydxTemp[6], yIn[7] ;

   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<nvar;i++) yIn[i]=yInput[i];

   G4double h = hstep * 0.5; 

   // Do two half steps

   DumbStepper(yIn,dydx,h,yTemp);
   RightHandSide(yTemp,dydxTemp) ;    
   DumbStepper(yTemp,dydxTemp,h,yOut); 

   // Store midpoint, chord calculation

   yMidPoint = G4ThreeVector( yTemp[0],  yTemp[1],  yTemp[2]); 

   // Do a full Step
   h = hstep ;
   DumbStepper(yIn,dydx,h,yTemp); 
   for(i=0;i<nvar;i++) {
      yErr[i] = yOut[i] - yTemp[i] ;
      yOut[i] += yErr[i]*correction ;    // Provides by 1 increased
                                         // order of accuracy
                                         // Richardson Extrapolation  
   }

   yInitial = G4ThreeVector( yIn[0],   yIn[1],   yIn[2]); 
   yFinal   = G4ThreeVector( yOut[0],  yOut[1],  yOut[2]); 

   return ;
}



G4double
G4MagErrorStepper::DistChord()   const 
{
  // Soon: must check whether h/R > 2 pi  !!
  //  Method below is good only for < 2 pi

  return G4LineSection::Distline( yMidPoint, yInitial, yFinal );
  // This is a class method that gives distance of Mid 
  //  from the Chord between the Initial and Final points.
}
 
