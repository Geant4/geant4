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
//  Tsitouras - 5(4) RK steppers   ( non-FSAL version )
//
//  Implements RK tableau from 'Table 1' of 
//    C. Tsitouras, “Runge–Kutta pairs of order 5(4) satisfying only
//    the first column simplifying assumption,”
//   Computers & Mathematics with Applications,
//   vol. 62, no. 2, pp. 770–775, 2011.
//
//  Adaptation / Geant4 implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 12 June 2015
//
// -------------------------------------------------------------------

#include "G4TsitourasRK45.hh"
#include "G4LineSection.hh"

/////////////////////////////////////////////////////////////////////
//
// Constructor

G4TsitourasRK45::G4TsitourasRK45(G4EquationOfMotion *EqRhs, 
				 G4int noIntegrationVariables, 
				 G4bool primary)
  : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
    fLastStepLength(0.), fAuxStepper(0)
{
  const G4int numberOfVariables = noIntegrationVariables;
  // G4cout << "G4TsitourasRK45 constructor called." << G4endl;
  
  ak2 = new G4double[numberOfVariables] ;  
  ak3 = new G4double[numberOfVariables] ; 
  ak4 = new G4double[numberOfVariables] ; 
  ak5 = new G4double[numberOfVariables] ; 
  ak6 = new G4double[numberOfVariables] ; 
  ak7 = new G4double[numberOfVariables] ;
  ak8 = new G4double[numberOfVariables] ;


  // Must ensure space extra 'state' variables exists - i.e. yIn[7]
  const G4int numStateMax  = std::max(GetNumberOfStateVariables(), 8);  
  const G4int numStateVars = std::max(noIntegrationVariables,
                                      numStateMax );
                                   // GetNumberOfStateVariables() ); 
                                      
  yTemp = new G4double[numStateVars] ;
  yIn = new G4double[numStateVars] ;

  fLastInitialVector = new G4double[numberOfVariables] ;
  fLastFinalVector = new G4double[numberOfVariables] ;

  fLastDyDx = new G4double[numberOfVariables];

  fMidVector = new G4double[numberOfVariables];
  fMidError =  new G4double[numberOfVariables];
  if( primary )
  { 
    fAuxStepper = new G4TsitourasRK45(EqRhs, numberOfVariables, !primary);
  }
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4TsitourasRK45::~G4TsitourasRK45()
{
  delete[] ak2;
  delete[] ak3;
  delete[] ak4;
  delete[] ak5;
  delete[] ak6;
  delete[] ak7;
  delete[] ak8;

  delete[] yTemp;
  delete[] yIn;

  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fLastDyDx;
  delete[] fMidVector;
  delete[] fMidError; 

  delete fAuxStepper;
}

//The following coefficients have been obtained from
// Table 1: The Coefficients of the new pair
//---Ref---
// C. Tsitouras, “Runge–Kutta pairs of order 5(4) satisfying only
// the first column simplifying assumption,”
// Computers & Mathematics with Applications,
// vol. 62, no. 2, pp. 770–775, 2011.
//-----------------------------------
// A corresponding matlab code was also found @ http://users.ntua.gr/tsitoura/new54.m

// Doing a step
void
G4TsitourasRK45::Stepper(  const G4double yInput[],
                         const G4double dydx[],
                               G4double Step,
                               G4double yOut[],
                               G4double yErr[])
{
    G4int i;

    const G4double
    b21 = 0.161 ,
    
    b31 = -0.00848065549235698854 ,
    b32 = 0.335480655492356989 ,
    
    b41 =  2.89715305710549343 ,
    b42 = -6.35944848997507484 ,
    b43 = 4.36229543286958141 ,

    b51 = 5.325864828439257,
    b52 = -11.748883564062828,
    b53 = 7.49553934288983621 ,
    b54 = -0.09249506636175525,

    b61 = 5.8614554429464200,
    b62 = -12.9209693178471093 ,
    b63 = 8.1593678985761586 ,
    b64 = -0.071584973281400997,
    b65 = -0.0282690503940683829,

    
    b71 = 0.0964607668180652295 ,
    b72 = 0.01, 
    b73 = 0.479889650414499575,
    b74 = 1.37900857410374189,
    b75 = -3.2900695154360807,
    b76 = 2.32471052409977398,
    
//    c1 = 0.001780011052226 ,
//    c2 = 0.000816434459657 ,
//    c3 = -0.007880878010262 ,
//    c4 = 0.144711007173263 ,
//    c5 = -0.582357165452555 ,
//    c6 = 0.458082105929187 ,
//    c7 = 1.0/66.0 ;

    dc1 = 0.0935237485818927066 - b71 , // - 0.001780011052226,
    dc2 = 0.00865288314156636761 - b72, // - 0.000816434459657,
    dc3 = 0.492893099131431868 - b73 ,  // + 0.007880878010262,
    dc4 = 1.14023541226785810 - b74 ,   //   0.144711007173263,
    dc5 = - 2.3291801924393646 - b75,   // + 0.582357165452555,
    dc6 = 1.56887504931661552 - b76 ,   // - 0.458082105929187,
    dc7 = 0.025; //- 1.0/66.0 ;
    
//    dc1 = -3.0/1280.0,
//    dc2 = 0.0,
//    dc3 = 6561.0/632320.0,
//    dc4 = -343.0/20800.0,
//    dc5 = 243.0/12800.0,
//    dc6 = -1.0/95.0,
//    dc7 = 0.0   ;
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];
    //  Saving yInput because yInput and yOut can be aliases for same array
    
    for(i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }

    // RightHandSide(yIn, dydx) ;
    // 1st Step - Not doing, getting passed
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
    }
    RightHandSide(yTemp, ak2) ;              // 2nd Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
    }
    RightHandSide(yTemp, ak3) ;              // 3rd Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    RightHandSide(yTemp, ak4) ;              // 4th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                                  b54*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5) ;              // 5th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b71*dydx[i] + b72*ak2[i] + b73*ak3[i] +
                                 b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yOut, ak7);       //7th Stage

    //Calculate the error in the step:
    for(i=0;i<numberOfVariables;i++)
    {
        yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i]  + dc6*ak6[i] + dc7*ak7[i] ) ;
        
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = dydx[i];
    }
    
    fLastStepLength = Step;
    
    return ;
}

void G4TsitourasRK45::SetupInterpolation() // (const G4double *yInput, const G4double *dydx, const G4double Step)
{
    //Nothing to be done
}


void G4TsitourasRK45::Interpolate(const G4double *yInput, const G4double *dydx, const G4double Step, G4double *yOut, G4double tau){

    
    G4double bf1, bf2, bf3, bf4, bf5, bf6, bf7;
    // Coefficients for all the seven stages.
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    G4double tau0 = tau;
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    G4double
    tau_2 = tau0*tau0 ;
//    tau_3 = tau0*tau_2,
//    tau_4 = tau_2*tau_2;

    bf1 = -1.0530884977290216*tau*(tau - 1.3299890189751412)*(tau_2 -
                1.4364028541716351*tau + 0.7139816917074209),
    bf2 = 0.1017*tau_2*(tau_2 - 2.1966568338249754*tau +
                        1.2949852507374631),
    bf3 = 2.490627285651252793*tau_2*(tau_2 - 2.38535645472061657*tau
                                      + 1.57803468208092486) ,
    bf4 = -16.54810288924490272*(tau - 1.21712927295533244)*
    				(tau - 0.61620406037800089)*tau_2,
    bf5 = 47.37952196281928122*(tau - 1.203071208372362603)*
    				(tau - 0.658047292653547382)*tau_2,
    bf6 = -34.87065786149660974*(tau - 1.2)*(tau -
                                0.666666666666666667)*tau_2,
    bf7 = 2.5*(tau - 1.0)*(tau - 0.6)*tau_2;
    
    //Putting together the coefficients calculated as the respective stage coefficients
    for( int i=0; i<numberOfVariables; i++){
        yOut[i] = yIn[i] + Step*(  bf1*dydx[i] + bf2*ak2[i] + bf3*ak3[i] + bf4*ak4[i]
                                 + bf5*ak5[i]  + bf6*ak6[i] + bf7*ak7[i]  ) ;
    }
}

///////////////////////////////////////////////////////////////////////////////

G4double  G4TsitourasRK45::DistChord() const
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
