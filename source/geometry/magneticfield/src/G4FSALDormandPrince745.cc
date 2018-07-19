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
//  DormandPrince7 - 5(4) implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 25 May 2015
//
//  G4FSALDormandPrince745.cc
//  Geant4
//
//    This is the source file of G4FSALDormandPrince745 class containing the
//    definition of the stepper() method that evaluates one step in
//    field propagation.
//    The Butcher table of the FDormand-Prince-7-4-5 method is as follows :
//
//    0   |
//    1/5 | 1/5
//    3/10| 3/40        9/40
//    4/5 | 44/45      −56/15      32/9
//    8/9 | 19372/6561 −25360/2187 64448/6561 −212/729
//    1   | 9017/3168  −355/33    46732/5247  49/176  −5103/18656
//    1   | 35/384      0         500/1113    125/192 −2187/6784    11/84
//    ---------------------------------------------------------------------------
//          35/384       0        500/1113    125/192  −2187/6784    11/84   0
//          5179/57600   0       7571/16695  393/640  −92097/339200 187/2100 1/40
//
//    Implementation by Somnath Banerjee - GSoC 2015
//       Work supported by Google as part of Google Summer of Code 2015.
//    Supervision / code review: John Apostolakis
//
//  First version: June 2015 - Somnath Banerjee

#include "G4FSALDormandPrince745.hh"
#include "G4LineSection.hh"
#include <cmath>

//Constructor
G4FSALDormandPrince745::G4FSALDormandPrince745(G4EquationOfMotion *EqRhs,
                                   G4int noIntegrationVariables,
                                   G4bool primary)
   : G4VFSALIntegrationStepper(EqRhs, noIntegrationVariables)
{
    
    const G4int numberOfVariables = noIntegrationVariables;
    
    //New Chunk of memory being created for use by the stepper
    
    //aki - for storing intermediate RHS
    ak2 = new G4double[numberOfVariables];
    ak3 = new G4double[numberOfVariables];
    ak4 = new G4double[numberOfVariables];
    ak5 = new G4double[numberOfVariables];
    ak6 = new G4double[numberOfVariables];
    ak7 = new G4double[numberOfVariables];
    // Also always allocate arrays for interpolation stages    
    ak8 = new G4double[numberOfVariables];
    ak9 = new G4double[numberOfVariables];
    
    yTemp = new G4double[numberOfVariables] ;
    yIn = new G4double[numberOfVariables] ;
    
    pseudoDydx_for_DistChord = new G4double[numberOfVariables];

    fInitialDyDx = new G4double[numberOfVariables];    
    fLastInitialVector = new G4double[numberOfVariables] ;
    fLastFinalVector = new G4double[numberOfVariables] ;
    fLastDyDx = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    
    fAuxStepper = nullptr;
    if( primary )
    {
        fAuxStepper = new G4FSALDormandPrince745(EqRhs, numberOfVariables,
                                           !primary);
    }
    fLastStepLength = -1.0;
}


//Destructor
G4FSALDormandPrince745::~G4FSALDormandPrince745()
{
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;  ak2=nullptr;
    delete[] ak3;  ak3=nullptr;
    delete[] ak4;  ak4=nullptr;
    delete[] ak5;  ak5=nullptr;
    delete[] ak6;  ak6=nullptr;
    delete[] ak7;  ak7=nullptr;
    delete[] ak8;  ak8=nullptr;
    delete[] ak9;  ak9=nullptr;
    
    delete[] yTemp; yTemp= nullptr;
    delete[] yIn;   yIn= nullptr;

    delete[] pseudoDydx_for_DistChord;  pseudoDydx_for_DistChord= nullptr;
    delete[] fInitialDyDx;              fInitialDyDx=       nullptr;
    
    delete[] fLastInitialVector;    fLastInitialVector= nullptr;
    delete[] fLastFinalVector;      fLastFinalVector  = nullptr;
    delete[] fLastDyDx;             fLastDyDx  = nullptr;
    delete[] fMidVector;            fMidVector = nullptr;
    delete[] fMidError;             fMidError  = nullptr;
    
    delete fAuxStepper;             fAuxStepper= nullptr;
}


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void G4FSALDormandPrince745::Stepper(const G4double yInput[],
                               const G4double dydx[],
                               G4double Step,
                               G4double yOut[],
                               G4double yErr[],
                               G4double nextDydx[]
                               )
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    b21 = 0.2 ,
    
    b31 = 3.0/40.0, b32 = 9.0/40.0 ,
    
    b41 = 44.0/45.0, b42 = -56.0/15.0, b43 = 32.0/9.0,
    
    b51 = 19372.0/6561.0, b52 = -25360.0/2187.0, b53 = 64448.0/6561.0,
    b54 = -212.0/729.0 ,
    
    b61 = 9017.0/3168.0 , b62 =   -355.0/33.0,
    b63 =  46732.0/5247.0    , b64 = 49.0/176.0 ,
    b65 = -5103.0/18656.0 ,
    
    b71 = 35.0/384.0, b72 = 0.,
    b73 = 500.0/1113.0, b74 = 125.0/192.0,
    b75 = -2187.0/6784.0, b76 = 11.0/84.0,
    
//    c1 = 35.0/384.0, c2 = .0,
//    c3 = 500.0/1113.0, c4 = 125.0/192.0,
//    c5 = -2187.0/6784.0, c6 = 11.0/84.0,
//    c7 = 0,
    
    dc1 = b71 - 5179.0/57600.0,
    dc2 = b72 - .0,
    dc3 = b73 - 7571.0/16695.0,
    dc4 = b74 - 393.0/640.0,
    dc5 = b75 + 92097.0/339200.0,
    dc6 = b76 - 187.0/2100.0,
    dc7 = - 1.0/40.0 ; //end of declaration
    
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    // The number of variables to be integrated over

    //  Saving yInput because yInput and yOut can be aliases for same array
    for(i=0;i<numberOfVariables;i++)
    {
        yIn[i]          = yInput[i];
        fInitialDyDx[i] = dydx[i];
    }
    // Ensure that time is initialised - in case it is not integrated
    yOut[7] = yTemp[7]  = yInput[7];
    
    // RightHandSide(yIn, DyDx) ;
    // 1st Step - Not doing, getting passed
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + b21*Step*fInitialDyDx[i] ;
    }
    RightHandSide(yTemp, ak2) ;              // 2nd Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b31*fInitialDyDx[i] + b32*ak2[i]) ;
    }
    RightHandSide(yTemp, ak3) ;              // 3rd Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b41*fInitialDyDx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    RightHandSide(yTemp, ak4) ;              // 4th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b51*fInitialDyDx[i] + b52*ak2[i] + b53*ak3[i] +
                                  b54*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5) ;              // 5th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*fInitialDyDx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i]          + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b71*fInitialDyDx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i]         + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yOut, ak7);               //7th and Final step
    
    for(i=0;i<numberOfVariables;i++)
    {
        
        yErr[i] = Step*(dc1*fInitialDyDx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i]          + dc6*ak6[i] + dc7*ak7[i] ) ;
        

        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = fInitialDyDx[i];
        nextDydx[i] = ak7[i];
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  G4FSALDormandPrince745::DistChord() const
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
                         fMidVector,   fMidError, pseudoDydx_for_DistChord );
    
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


void G4FSALDormandPrince745::interpolate(  const G4double yInput[],
                                     const G4double dydx[],
                                     G4double yOut[],
                                     G4double Step,
                                      G4double tau){
    
    G4double
    bf1, bf2, bf3, bf4, bf5, bf6, bf7;

    
    

    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    G4double tau0 = tau;
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    G4double
    tau_2 = tau0*tau0 ,
    tau_3 = tau0*tau_2,
    tau_4 = tau_2*tau_2;
    
    bf1 = (157015080.0*tau_4 - 13107642775.0*tau_3+ 34969693132.0*tau_2- 32272833064.0*tau
           + 11282082432.0)/11282082432.0,
    bf2 = 0.0 ,
    bf3 = - 100.0*tau*(15701508.0*tau_3 - 914128567.0*tau_2 + 2074956840.0*tau
                 - 1323431896.0)/32700410799.0,
    bf4 = 25.0*tau*(94209048.0*tau_3- 1518414297.0*tau_2+ 2460397220.0*tau - 889289856.0)/5641041216.0 ,
    bf5 = -2187.0*tau*(52338360.0*tau_3 - 451824525.0*tau_2 + 687873124.0*tau - 259006536.0)/199316789632.0 ,
    bf6 =  11.0*tau*(106151040.0*tau_3- 661884105.0*tau_2 + 946554244.0*tau - 361440756.0)/2467955532.0 ,
    bf7 = tau*(1.0 - tau)*(8293050.0*tau_2 - 82437520.0*tau + 44764047.0)/ 29380423.0 ;

    
    for( int i=0; i<numberOfVariables; i++){
        yOut[i] = yIn[i] + Step*tau*(bf1*dydx[i] + bf2*ak2[i] + bf3*ak3[i] + bf4*ak4[i]
                                     + bf5*ak5[i] + bf6*ak6[i] + bf7*ak7[i]  ) ;
    }
    

    
}

void G4FSALDormandPrince745::SetupInterpolate(const G4double yInput[],
                                   const G4double dydx[],
                                   const G4double Step ){
    
    //Coefficients for the additional stages :
    G4double
    b81 =  6245.0/62208.0 ,
    b82 =  0.0 ,
    b83 =  8875.0/103032.0 ,
    b84 = -125.0/1728.0 ,
    b85 =  801.0/13568.0 ,
    b86 = -13519.0/368064.0 ,
    b87 =  11105.0/368064.0 ,
    
    b91 =  632855.0/4478976.0 ,
    b92 =  0.0 ,
    b93 =  4146875.0/6491016.0 ,
    b94 =  5490625.0/14183424.0 ,
    b95 = -15975.0/108544.0 ,
    b96 =  8295925.0/220286304.0 ,
    b97 = -1779595.0/62938944.0 ,
    b98 = -805.0/4104.0 ;
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    //  Saving yInput because yInput and yOut can be aliases for same array
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    yTemp[7]  = yIn[7];
    
    //Evaluate the extra stages :
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*( b81*dydx[i] + b82*ak2[i] + b83*ak3[i] +
                                   b84*ak4[i] + b85*ak5[i] + b86*ak6[i] +
                                   b87*ak7[i] );
    }
    RightHandSide( yTemp, ak8 );              //8th Stage
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step * ( b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                   b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
                                   b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide( yTemp, ak9 );          //9th Stage
    
    
    
}



void G4FSALDormandPrince745::Interpolate( const G4double yInput[],
                 const G4double dydx[],
                 const G4double Step,
                 G4double yOut[],
                                G4double tau ){
    //Define the coefficients for the polynomials
    G4double bi[10][5], b[10];
    G4int numberOfVariables = this->GetNumberOfVariables();
    
    //  COEFFICIENTS OF   bi[1]
    bi[1][0] =  1.0 ,
    bi[1][1] = -38039.0/7040.0 ,
    bi[1][2] =  125923.0/10560.0 ,
    bi[1][3] = -19683.0/1760.0 ,
    bi[1][4] =  3303.0/880.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[2]
    bi[2][0] =  0.0 ,
    bi[2][1] =  0.0 ,
    bi[2][2] =  0.0 ,
    bi[2][3] =  0.0 ,
    bi[2][4] =  0.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[3]
    bi[3][0] =  0.0 ,
    bi[3][1] = -12500.0/4081.0 ,
    bi[3][2] =  205000.0/12243.0 ,
    bi[3][3] = -90000.0/4081.0 ,
    bi[3][4] =  36000.0/4081.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[4]
    bi[4][0] =  0.0 ,
    bi[4][1] = -3125.0/704.0 ,
    bi[4][2] =  25625.0/1056.0 ,
    bi[4][3] = -5625.0/176.0 ,
    bi[4][4] =  1125.0/88.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[5]
    bi[5][0] =  0.0 ,
    bi[5][1] =  164025.0/74624.0 ,
    bi[5][2] = -448335.0/37312.0 ,
    bi[5][3] =  295245.0/18656.0 ,
    bi[5][4] = -59049.0/9328.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[6]
    bi[6][0] =  0.0 ,
    bi[6][1] = -25.0/28.0 ,
    bi[6][2] =  205.0/42.0 ,
    bi[6][3] = -45.0/7.0 ,
    bi[6][4] =  18.0/7.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[7]
    bi[7][0] =  0.0 ,
    bi[7][1] = -2.0/11.0 ,
    bi[7][2] =  73.0/55.0 ,
    bi[7][3] = -171.0/55.0 ,
    bi[7][4] =  108.0/55.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[8]
    bi[8][0] =  0.0 ,
    bi[8][1] =  189.0/22.0 ,
    bi[8][2] = -1593.0/55.0 ,
    bi[8][3] =  3537.0/110.0 ,
    bi[8][4] = -648.0/55.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[9]
    bi[9][0] =  0.0 ,
    bi[9][1] =  351.0/110.0 ,
    bi[9][2] = -999.0/55.0 ,
    bi[9][3] =  2943.0/110.0 ,
    bi[9][4] = -648.0/55.0 ;
    //  --------------------------------------------------------
    

    
    for(G4int i = 0; i< numberOfVariables; i++)
        yIn[i] = yInput[i];
    
    G4double tau0 = tau;
    //    Calculating the polynomials :
    
    for(int i=1; i<=9; i++){    //Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = 1.0;
        for(int j=0; j<=4; j++){
            b[i] += bi[i][j]*tau;
            tau*=tau0;
        }
    }
    
    for(int i=0; i<numberOfVariables; i++){     //Here i IS the cooridnate no.
        yOut[i] = yIn[i] + Step*tau0*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                      b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                      b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] );
    }

}







