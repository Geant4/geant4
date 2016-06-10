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
// $Id: G4ConstRK4.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// - 18.09.2008 - J.Apostolakis, T.Nikitina - Created
// -------------------------------------------------------------------

#include "G4ConstRK4.hh"
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

//////////////////////////////////////////////////////////////////
//
// Constructor sets the number of *State* variables (default = 8)
//   The number of variables integrated is always 6

G4ConstRK4::G4ConstRK4(G4Mag_EqRhs* EqRhs, G4int numStateVariables)
  : G4MagErrorStepper(EqRhs, 6, numStateVariables)
{
  // const G4int numberOfVariables= 6;
  if( numStateVariables < 8 ) 
  {
    std::ostringstream message;
    message << "The number of State variables at least 8 " << G4endl
            << "Instead it is - numStateVariables= " << numStateVariables;
    G4Exception("G4ConstRK4::G4ConstRK4()", "GeomField0002",
                FatalException, message, "Use another Stepper!");
  }

  fEq = EqRhs;
  yMiddle  = new G4double[8];
  dydxMid  = new G4double[8];
  yInitial = new G4double[8];
  yOneStep = new G4double[8];

  dydxm = new G4double[8];
  dydxt = new G4double[8]; 
  yt    = new G4double[8]; 
  Field[0]=0.; Field[1]=0.; Field[2]=0.;
}

////////////////////////////////////////////////////////////////
//
// Destructor

G4ConstRK4::~G4ConstRK4()
{
   delete [] yMiddle;
   delete [] dydxMid;
   delete [] yInitial;
   delete [] yOneStep;
   delete [] dydxm;
   delete [] dydxt;
   delete [] yt;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for the variables y[0,..,n-1] and their derivatives
// dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
// method to advance the solution over an interval h and return the
// incremented variables as yout[0,...,n-1], which is not a distinct
// array from y. The user supplies the routine RightHandSide(x,y,dydx),
// which returns derivatives dydx at x. The source is routine rk4 from
// NRC p. 712-713 .

void G4ConstRK4::DumbStepper( const G4double  yIn[],
                              const G4double  dydx[],
                                    G4double  h,
                                    G4double  yOut[])
{
   G4double  hh = h*0.5 , h6 = h/6.0  ;
   
   // 1st Step K1=h*dydx
   yt[5] = yIn[5] + hh*dydx[5] ;
   yt[4] = yIn[4] + hh*dydx[4] ;
   yt[3] = yIn[3] + hh*dydx[3] ;
   yt[2] = yIn[2] + hh*dydx[2] ;
   yt[1] = yIn[1] + hh*dydx[1] ;
   yt[0] = yIn[0] + hh*dydx[0] ;
   RightHandSideConst(yt,dydxt) ;        

   // 2nd Step K2=h*dydxt
   yt[5] = yIn[5] + hh*dydxt[5] ;
   yt[4] = yIn[4] + hh*dydxt[4] ;
   yt[3] = yIn[3] + hh*dydxt[3] ;
   yt[2] = yIn[2] + hh*dydxt[2] ;
   yt[1] = yIn[1] + hh*dydxt[1] ;
   yt[0] = yIn[0] + hh*dydxt[0] ;
   RightHandSideConst(yt,dydxm) ;     

   // 3rd Step K3=h*dydxm
   // now dydxm=(K2+K3)/h
   yt[5] = yIn[5] + h*dydxm[5] ;
   dydxm[5] += dydxt[5] ;  
   yt[4] = yIn[4] + h*dydxm[4] ;
   dydxm[4] += dydxt[4] ;  
   yt[3] = yIn[3] + h*dydxm[3] ;
   dydxm[3] += dydxt[3] ;  
   yt[2] = yIn[2] + h*dydxm[2] ;
   dydxm[2] += dydxt[2] ;  
   yt[1] = yIn[1] + h*dydxm[1] ;
   dydxm[1] += dydxt[1] ;  
   yt[0] = yIn[0] + h*dydxm[0] ;
   dydxm[0] += dydxt[0] ;  
   RightHandSideConst(yt,dydxt) ;   

   // 4th Step K4=h*dydxt
   yOut[5] = yIn[5]+h6*(dydx[5]+dydxt[5]+2.0*dydxm[5]);
   yOut[4] = yIn[4]+h6*(dydx[4]+dydxt[4]+2.0*dydxm[4]);
   yOut[3] = yIn[3]+h6*(dydx[3]+dydxt[3]+2.0*dydxm[3]);
   yOut[2] = yIn[2]+h6*(dydx[2]+dydxt[2]+2.0*dydxm[2]);
   yOut[1] = yIn[1]+h6*(dydx[1]+dydxt[1]+2.0*dydxm[1]);
   yOut[0] = yIn[0]+h6*(dydx[0]+dydxt[0]+2.0*dydxm[0]);
   
}  // end of DumbStepper ....................................................

////////////////////////////////////////////////////////////////
//
// Stepper

void
G4ConstRK4::Stepper( const G4double yInput[],
                     const G4double dydx[],
                           G4double hstep,
                           G4double yOutput[],
                           G4double yError [] )
{
   const G4int nvar = 6;  // number of variables integrated
   const G4int maxvar= GetNumberOfStateVariables();

   // Correction for Richardson extrapolation
   G4double  correction = 1. / ( (1 << IntegratorOrder()) -1 );

   G4int i;
   
   // Saving yInput because yInput and yOutput can be aliases for same array
   for (i=0;    i<maxvar; i++) { yInitial[i]= yInput[i]; }
 
   // Must copy the part of the state *not* integrated to the output
   for (i=nvar; i<maxvar; i++) { yOutput[i]=  yInput[i]; }

   // yInitial[7]= yInput[7];  //  The time is typically needed
   yMiddle[7]  = yInput[7];   // Copy the time from initial value 
   yOneStep[7] = yInput[7];   // As it contributes to final value of yOutput ?
   // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
   yError[7] = 0.0;         

   G4double halfStep = hstep * 0.5; 

   // Do two half steps
   //
   GetConstField(yInitial,Field);
   DumbStepper  (yInitial,  dydx,   halfStep, yMiddle);
   RightHandSideConst(yMiddle, dydxMid);    
   DumbStepper  (yMiddle, dydxMid, halfStep, yOutput); 

   // Store midpoint, chord calculation
   //
   fMidPoint = G4ThreeVector( yMiddle[0],  yMiddle[1],  yMiddle[2]); 

   // Do a full Step
   //
   DumbStepper(yInitial, dydx, hstep, yOneStep);
   for(i=0;i<nvar;i++)
   {
      yError [i] = yOutput[i] - yOneStep[i] ;
      yOutput[i] += yError[i]*correction ;
        // Provides accuracy increased by 1 order via the 
        // Richardson extrapolation  
   }

   fInitialPoint = G4ThreeVector( yInitial[0], yInitial[1], yInitial[2]); 
   fFinalPoint   = G4ThreeVector( yOutput[0],  yOutput[1],  yOutput[2]); 

   return;
}

////////////////////////////////////////////////////////////////
//
// Estimate the maximum distance from the curve to the chord
//
// We estimate this using the distance of the midpoint to chord.
// The method below is good only for angle deviations < 2 pi;
// this restriction should not be a problem for the Runge Kutta methods, 
// which generally cannot integrate accurately for large angle deviations

G4double G4ConstRK4::DistChord() const 
{
  G4double distLine, distChord; 

  if (fInitialPoint != fFinalPoint)
  {
     distLine= G4LineSection::Distline( fMidPoint, fInitialPoint, fFinalPoint );
       // This is a class method that gives distance of Mid 
       // from the Chord between the Initial and Final points
     distChord = distLine;
  }
  else
  {
     distChord = (fMidPoint-fInitialPoint).mag();
  }
  return distChord;
}
