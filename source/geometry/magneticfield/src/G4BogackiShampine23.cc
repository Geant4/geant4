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
//  Bogacki-Shampine - 4 - 3(2) non-FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 20 May 2015
//
//  History
// -----------------------------
//  Created by Somnath Banerjee on 20 May 2015
///////////////////////////////////////////////////////////////////////////////


/*

This contains the stepper function of the G4BogackiShampine23 class

The Bogacki shampine method has the following Butcher's tableau

0  |
1/2|1/2
3/4|0	3/4
1  |2/9	1/3	4/9
-------------------
   |2/9	1/3	4/9	0
   |7/24 1/4 1/3 1/8

*/

#include "G4BogackiShampine23.hh"
#include "G4LineSection.hh"

// using namespace std;

//Constructor
G4BogackiShampine23::G4BogackiShampine23(G4EquationOfMotion *EqRhs,
				 G4int noIntegrationVariables,
				 G4bool primary)
  : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
    fLastStepLength(0.), fAuxStepper(nullptr)
{
  const G4int numberOfVariables = noIntegrationVariables;

  SetIntegrationOrder(2);  // Apparently 3-rd order extension not used
  SetFSAL(true);

  ak2 = new G4double[numberOfVariables] ;
  ak3 = new G4double[numberOfVariables] ;
  ak4 = new G4double[numberOfVariables] ;

  pseudoDydx_for_DistChord = new G4double[numberOfVariables];

  const G4int numStateVars = std::max(noIntegrationVariables,
                                      GetNumberOfStateVariables() );  

  yTemp = new G4double[numberOfVariables] ;
  yIn = new G4double[numberOfVariables] ;
  
  fLastInitialVector = new G4double[numStateVars] ;
  fLastFinalVector = new G4double[numStateVars] ;
  fLastDyDx = new G4double[numStateVars];

  fMidVector = new G4double[numStateVars];
  fMidError =  new G4double[numStateVars];
  if( primary )
  {
    fAuxStepper = new G4BogackiShampine23(EqRhs, numberOfVariables, !primary);
  }
}


//Destructor
G4BogackiShampine23::~G4BogackiShampine23()
{
  delete[] ak2;
  delete[] ak3;
  delete[] ak4;
  delete[] pseudoDydx_for_DistChord;

  delete[] yTemp;
  delete[] yIn;

  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fLastDyDx;
  delete[] fMidVector;
  delete[] fMidError;

  delete fAuxStepper;
}

//******************************************************************************
//
// Given values for n = 4 variables yIn[0,...,n-1]
// known  at x, use the 3rd order Bogacki Shampine method
// to advance the solution over an interval Step
// and return the incremented variables as yOut[0,...,n-1]. Also
// return an estimate of the local truncation error yErr[] using the
// embedded 2nd order method. The user supplies routine
// RightHandSide(y,dydx), which returns derivatives dydx for y .


//******************************************************************************


void
G4BogackiShampine23::Stepper( const G4double yInput[],
                         	const G4double DyDx[],
                              	  G4double Step,
                              	  G4double yOut[],
                              	  G4double yErr[])
{
 G4int i;

 const G4double  b21 = 0.5 ,
                 b31 = 0. , b32 = 3.0/4.0 ,
                 b41 = 2.0/9.0, b42 = 1.0/3.0 , b43 = 4.0/9.0;

 const G4double  dc1 = b41 - 7.0/24.0 ,  dc2 = b42 - 1.0/4.0 ,
  				 dc3 = b43 - 1.0/3.0 , dc4 = - 0.125 ;
    
 // Initialise time to t0, needed when it is not updated by the integration.
 //        [ Note: Only for time dependent fields (usually electric)
 //                  is it neccessary to integrate the time.]
 yOut[7] = yTemp[7]   = yIn[7];

 const G4int numberOfVariables= this->GetNumberOfVariables();	// The number of variables to be integrated over

   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<numberOfVariables;i++)
   {
     	yIn[i]=yInput[i];
   }
    // RightHandSide(yIn, dydx) ;              // 1st Step --Not doing, getting passed
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + b21*Step*DyDx[i] ;
    }
    RightHandSide(yTemp, ak2) ;              // 2nd Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b31*DyDx[i] + b32*ak2[i]) ;
    }
    RightHandSide(yTemp, ak3) ;              // 3rd Step

    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b41*DyDx[i] + b42*ak2[i] + b43*ak3[i]) ;
        // yOut[i] = yIn[i] + Step*(c1*DyDx[i]+ c2*ak2[i] + c3*ak3[i] + c4*ak4[i]);
    }
    // Extra step used only in calculation of error
    RightHandSide(yOut, ak4) ;              // 4th Step
    // Derivative and end-point already calculated in 'ak4' ! => Can be used in FSAL version
    
    for(i=0;i<numberOfVariables;i++)
    {
        yErr[i] = Step*(dc1*DyDx[i] + dc2*ak2[i] + dc3*ak3[i] +
                        dc4*ak4[i] ) ;
        
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = DyDx[i];
    }
    // NormaliseTangentVector( yOut ); // Not wanted
    
    fLastStepLength =Step;
    
    return ;
}

G4double  G4BogackiShampine23::DistChord() const
{
  G4double distLine, distChord;
  G4ThreeVector initialPoint, finalPoint, midPoint;

  // Store last initial and final points (they will be overwritten in self-Stepper call!)
  initialPoint = G4ThreeVector( fLastInitialVector[0],
                                fLastInitialVector[1], fLastInitialVector[2]);
  finalPoint   = G4ThreeVector( fLastFinalVector[0],
                                fLastFinalVector[1],  fLastFinalVector[2]);

  // Do half a step using StepNoErr

  fAuxStepper->Stepper( fLastInitialVector, fLastDyDx, 0.5 * fLastStepLength,
           fMidVector,   fMidError );

  midPoint = G4ThreeVector( fMidVector[0], fMidVector[1], fMidVector[2]);

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord


  if (initialPoint != finalPoint)
  {
     distLine  = G4LineSection::Distline( midPoint, initialPoint, finalPoint );
     distChord = distLine;
  }
  else
  {
     distChord = (midPoint-initialPoint).mag();
  }
  return distChord;
}

//------Verified-------
