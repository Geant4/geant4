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
//  Dormand-Lockyer-McGorrigan-Prince-6-3-4 non-FSAL implementation
//     RK4(3)6FD - forced Non-FSAL
//
//  Design/implementation  by Somnath Banerjee
//     Sponsored by Google in Google Summer of Code 2015.
//  Supervision / code review: John Apostolakis
//
// First version: 7 July 2015
//
//  G4DoLoMcPriRK34.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 7 July 2015
//
//    This is the source file of G4DoLoMcPriRK34 class containing the
//    definition of the Stepper() method that evaluates one Step in
//    field propagation.
//
//    The Butcher table of the Dormand-Lockyer-McGorrigan-Prince-6-3-4 method is as follows :
//    [ to be added here ]

#include "G4DoLoMcPriRK34.hh"
#include "G4LineSection.hh"

// Constructor
G4DoLoMcPriRK34::G4DoLoMcPriRK34(G4EquationOfMotion *EqRhs,
                                  G4int noIntegrationVariables,
                                  G4bool primary)
   : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
     fLastStepLength( -1.0 ), fAuxStepper( nullptr )     
{
   const G4int numberOfVariables = noIntegrationVariables;
   
   //New Chunk of memory being created for use by the Stepper
   
   //aki - for storing intermediate RHS
   ak2 = new G4double[numberOfVariables];
   ak3 = new G4double[numberOfVariables];
   ak4 = new G4double[numberOfVariables];
   ak5 = new G4double[numberOfVariables];
   ak6 = new G4double[numberOfVariables];
   
   yTemp = new G4double[numberOfVariables] ;
   yIn = new G4double[numberOfVariables] ;
      
   fLastInitialVector = new G4double[numberOfVariables] ;
   fLastFinalVector = new G4double[numberOfVariables] ;
   fLastDyDx = new G4double[numberOfVariables];
   
   fMidVector = new G4double[numberOfVariables];
   fMidError =  new G4double[numberOfVariables];
   if( primary )
   {
       fAuxStepper = new G4DoLoMcPriRK34(EqRhs, numberOfVariables,
                                          !primary);
   }
}


//Destructor
G4DoLoMcPriRK34::~G4DoLoMcPriRK34()
{
   //clear all previously allocated memory for Stepper and DistChord
   delete[] ak2;
   delete[] ak3;
   delete[] ak4;
   delete[] ak5;
   delete[] ak6;
   
   delete[] yTemp;
   delete[] yIn;
   
   delete[] fLastInitialVector;
   delete[] fLastFinalVector;
   delete[] fLastDyDx;
   delete[] fMidVector;
   delete[] fMidError;
   
   delete fAuxStepper;
   
   
}


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void G4DoLoMcPriRK34::Stepper(const G4double yInput[],
                          const G4double DyDx[],
                          	    G4double Step,
                            	G4double yOut[],
                              	G4double yErr[] )
{
   G4int i;
   
   //The various constants defined on the basis of butcher tableu
   const G4double  //G4double - only once
   
   b21 = 7.0/27.0 , 
   

   b31 = 7.0/72.0 , 
   b32 = 7.0/24.0 ,
   
   b41 = 3043.0/3528.0 , 
   b42 = -3757.0/1176.0 ,
   b43 = 1445.0/441.0,
   
   b51 =  17617.0/11662.0 ,
   b52 = -4023.0/686.0 , 
   b53 =  9372.0/1715.0 ,
   b54 = -66.0/595.0 ,
   
   b61 =  29.0/238.0 ,
   b62 =  0.0 , 
   b63 =  216.0/385.0 ,
   b64 =  54.0/85.0 ,
   b65 =  -7.0/22.0 ,
   

   
   dc1 = 363.0/2975.0 - b61 ,
   dc2 = 0.0 - b62 ,
   dc3 = 981.0/1750.0 - b63,
   dc4 = 2709.0/4250.0 - b64 ,
   dc5 = -3.0/10.0 - b65 ,
   dc6 = -1.0/50.0 ; //end of declaration
   
   
   const G4int numberOfVariables= this->GetNumberOfVariables();
   
   // The number of variables to be integrated over
   yOut[7] = yTemp[7]  = yIn[7];
   //  Saving yInput because yInput and yOut can be aliases for same array
   
   for(i=0;i<numberOfVariables;i++)
   {
       yIn[i]=yInput[i];
   }
   
   
   
   // RightHandSide(yIn, DyDx) ;
   // 1st stage - Not doing, getting passed
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + b21*Step*DyDx[i] ;
   }
   RightHandSide(yTemp, ak2) ;              // 2nd stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + Step*(b31*DyDx[i] + b32*ak2[i]) ;
   }
   RightHandSide(yTemp, ak3) ;              // 3rd stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + Step*(b41*DyDx[i] + b42*ak2[i] + b43*ak3[i]) ;
   }
   RightHandSide(yTemp, ak4) ;              // 4th stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + Step*(b51*DyDx[i] + b52*ak2[i] + b53*ak3[i] +
                                 b54*ak4[i]) ;
   }
   RightHandSide(yTemp, ak5) ;              // 5th stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yOut[i] = yIn[i] + Step*(b61*DyDx[i] + b62*ak2[i] + b63*ak3[i] +
                                 b64*ak4[i] + b65*ak5[i]) ;
   }
   RightHandSide(yOut, ak6) ;              // 6th and Final stage
   

   
   for(i=0;i<numberOfVariables;i++)
   {
       
       yErr[i] = Step*(dc1*DyDx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                       dc5*ak5[i] + dc6*ak6[i] ) ;
       

       // Store Input and Final values, for possible use in calculating chord
       fLastInitialVector[i] = yIn[i] ;
       fLastFinalVector[i]   = yOut[i];
       fLastDyDx[i]          = DyDx[i];
       
       
   }
   
   fLastStepLength = Step;
   
   return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  G4DoLoMcPriRK34::DistChord() const
{
   G4double distLine, distChord;
   G4ThreeVector initialPoint, finalPoint, midPoint;
   
   // Store last initial and final points (they will be overwritten in self-Stepper call!)
   initialPoint = G4ThreeVector( fLastInitialVector[0],
                                fLastInitialVector[1], fLastInitialVector[2]);
   finalPoint   = G4ThreeVector( fLastFinalVector[0],
                                fLastFinalVector[1],  fLastFinalVector[2]);
   
   // Do half a Step using StepNoErr
   
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

void G4DoLoMcPriRK34::SetupInterpolation()
{}

void G4DoLoMcPriRK34::SetupInterpolate( const G4double /* yInput */ [] ,
                                      const G4double  /* dydx */ [] ,
                                      const G4double  /* Step */ )
{
    //Do Nothing
}


void G4DoLoMcPriRK34::Interpolate( G4double tau, 
                                 G4double yOut[])
{
   Interpolate( fLastInitialVector, fLastDyDx, fLastStepLength, yOut, tau );
}

// Function to evaluate the interpolation at tau fraction of the step
void G4DoLoMcPriRK34::Interpolate( const G4double yInput[],
                                const G4double dydx[],
                                const G4double Step,
                                G4double yOut[],
                                G4double tau ){
   G4double
   bf1, bf2, bf3, bf4, bf5, bf6;    

   
   const G4int numberOfVariables= this->GetNumberOfVariables();
   
   for(int i=0;i<numberOfVariables;i++)
   {
       yIn[i]=yInput[i];
   }
   
   G4double
   	tau_2 = tau*tau ,
	tau_3 = tau*tau_2;
   
    //Calculating the polynomials (coefficients for the respective stages)
   bf1 = -(162.0*tau_3 - 504.0*tau_2 + 551.0*tau - 238.0)/238.0 ,
   bf2 =  0.0 ,
   bf3 =  27.0*tau*(27.0*tau_2 - 70.0*tau + 51.0 )/385.0 ,
   bf4 = -27*tau*(27.0*tau_2 - 50.0*tau + 21.0)/85.0 ,
   bf5 =  7.0*tau*(2232.0*tau_2 - 4166.0*tau + 1785.0 )/3278.0 ,
   bf6 = tau*(tau - 1.0)*(387.0*tau - 238.0)/149.0 ;
   
   for( int i=0; i<numberOfVariables; i++){
       yOut[i] = yIn[i] + Step*tau*(bf1*dydx[i] + bf2*ak2[i] + bf3*ak3[i] + 
								       	bf4*ak4[i] + bf5*ak5[i] + bf6*ak6[i] ) ;
   }
   

   
}


//-------Verified------- - hackabot

