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
// $Id: G4DormandPrince745.cc 107470 2017-11-15 07:14:28Z gcosmo $
//
// Class description:
//
//  DormandPrince7 - 5(4) non-FSAL
//
//    This is the source file of G4DormandPrince745 class containing the
//    definition of the stepper() method that evaluates one step in
//    field propagation.
//	  The coefficients and the algorithm have been adapted from
//
//    Table 2 : Coefficients of RK5(4)7M
//	  ---Ref---
//    J. R. Dormand and P. J. Prince, “A family of embedded Runge-Kutta formulae,”
//		Journal of computational and applied …, vol. 6, no. 1, pp. 19–26, 1980.
//		------------------
//
//    The Butcher table of the Dormand-Prince-7-4-5 method is as follows :
//
//    0   |
//    1/5 | 1/5
//    3/10| 3/40        9/40
//    4/5 | 44/45      −56/15      32/9
//    8/9 | 19372/6561 −25360/2187 64448/6561 −212/729
//    1   | 9017/3168  −355/33    46732/5247  49/176  −5103/18656
//    1   | 35/384      0         500/1113    125/192 −2187/6784    11/84
//    ------------------------------------------------------------------------
//          35/384       0        500/1113    125/192  −2187/6784    11/84   0
//          5179/57600   0       7571/16695  393/640  −92097/339200 187/2100 1/40
//
//
//    Implementation by Somnath Banerjee - GSoC 2015
//       Work supported by Google as part of Google Summer of Code 2015.
//    Supervision / code review: John Apostolakis
//
//  First version: 25 May 2015 - Somnath Banerjee
//
//  Note: Current version includes 3 versions of 'DistChord' method.
//        Default is hard-coded interpolation.
//
#include "G4DormandPrince745.hh"
#include "G4LineSection.hh"
#include <cmath>

//Constructor
G4DormandPrince745::G4DormandPrince745(G4EquationOfMotion *EqRhs,
                                   G4int noIntegrationVariables,
                                   G4bool primary)
: G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
  fAuxStepper(0)  
{
    const G4int numberOfVariables = // noIntegrationVariables;
      std::max( noIntegrationVariables,
               ( ( (noIntegrationVariables-1)/4 + 1 ) * 4 ) );
    // For better alignment with cache-line
    
    //New Chunk of memory being created for use by the stepper
    
    //ak_i - for storing intermediate RHS
    ak2 = new G4double[numberOfVariables];
    ak3 = new G4double[numberOfVariables];
    ak4 = new G4double[numberOfVariables];
    ak5 = new G4double[numberOfVariables];
    ak6 = new G4double[numberOfVariables];
    ak7 = new G4double[numberOfVariables];
    // Also always allocate arrays for interpolation stages
    ak8 = new G4double[numberOfVariables];
    ak9 = new G4double[numberOfVariables];

    // Must ensure space for extra 'state' variables exists - i.e. yIn[7]
    const G4int numStateVars =
       std::max(noIntegrationVariables,
                std::max( GetNumberOfStateVariables(), 8)
               );
    yTemp = new G4double[numStateVars];
    yIn   = new G4double[numStateVars];
    
    fLastInitialVector = new G4double[numStateVars] ;
    fLastFinalVector = new G4double[numStateVars] ;

    // fLastDyDx  = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numStateVars];
    fMidError =  new G4double[numStateVars];
    
    yTemp = new G4double[numberOfVariables] ;
    yIn =   new G4double[numberOfVariables] ;

    fLastInitialVector = new G4double[numberOfVariables] ;
    fLastFinalVector = new G4double[numberOfVariables] ;
    fInitialDyDx = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    fAuxStepper = nullptr;
    if( primary )
    {
        fAuxStepper = new G4DormandPrince745(EqRhs, numberOfVariables,
                                           !primary);
    }
    fLastStepLength = -1.0;
}

//Destructor
G4DormandPrince745::~G4DormandPrince745()
{
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;
    delete[] ak3;
    delete[] ak4;
    delete[] ak5;
    delete[] ak6;
    delete[] ak7;
    // Used only for interpolation
    delete[] ak8;
    delete[] ak9;
    
    delete[] yTemp;
    delete[] yIn;
    
    delete[] fLastInitialVector;
    delete[] fLastFinalVector;
    delete[] fInitialDyDx;
    delete[] fMidVector;
    delete[] fMidError;
    
    delete fAuxStepper;
}


//	  The coefficients and the algorithm have been adapted from
//    Table 2 : Coefficients of RK5(4)7M
//	  ---Ref---
//    J. R. Dormand and P. J. Prince, “A family of embedded Runge-Kutta formulae,”
//		Journal of computational and applied …, vol. 6, no. 1, pp. 19–26, 1980.
//		------------------

//    The Butcher table of the Dormand-Prince-7-4-5 method is as follows :
//
//    0   |
//    1/5 | 1/5
//    3/10| 3/40        9/40
//    4/5 | 44/45      −56/15      32/9
//    8/9 | 19372/6561 −25360/2187 64448/6561 −212/729
//    1   | 9017/3168  −355/33    46732/5247  49/176  −5103/18656
//    1   | 35/384      0         500/1113    125/192 −2187/6784    11/84
//    ------------------------------------------------------------------------
//          35/384       0        500/1113    125/192  −2187/6784    11/84   0
//          5179/57600   0       7571/16695  393/640  −92097/339200 187/2100 1/40


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void G4DormandPrince745::Stepper(const G4double yInput[],
                               const G4double DyDx[],
                                     G4double Step,
                                     G4double yOut[],
                                     G4double yErr[] )
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
    
    //Sum of columns, sum(bij) = ei
//    e1 = 0. ,
//    e2 = 1.0/5.0 ,
//    e3 = 3.0/10.0 ,
//    e4 = 4.0/5.0 ,
//    e5 = 8.0/9.0 ,
//    e6 = 1.0 ,
//    e7 = 1.0 ,
    
// Difference between the higher and the lower order method coeff. :
    // b7j are the coefficients of higher order
    
    dc1 = -( b71 - 5179.0/57600.0),
    dc2 = -( b72 - .0),
    dc3 = -( b73 - 7571.0/16695.0),
    dc4 = -( b74 - 393.0/640.0),
    dc5 = -( b75 + 92097.0/339200.0),
    dc6 = -( b76 - 187.0/2100.0),
    dc7 = -( - 1.0/40.0 ); //end of declaration
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yInput[7];
    //  Saving yInput because yInput and yOut can be aliases for same array
    
    for(i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    // RightHandSide(yIn, DyDx) ;
    // 1st Step - Not doing, getting passed
    
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
        yTemp[i] = yIn[i] + Step*(b61*DyDx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b71*DyDx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i] + b75*ak5[i] + b76*ak6[i] );
    }
    RightHandSide(yOut, ak7);				//7th and Final stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yErr[i] = Step*(dc1*DyDx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                            dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] ) + 1.5e-18 ;

        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fInitialDyDx[i]          = DyDx[i];        
    }
    
    fLastStepLength = Step;
    
    return ;
}


// Calculate DistChord given start, mid and end-point of step
G4double G4DormandPrince745::DistLine( G4double yStart[], G4double yMid[], G4double yEnd[] ) const
{
    G4double distLine, distChord;
    G4ThreeVector initialPoint, finalPoint, midPoint;
    
    initialPoint = G4ThreeVector( yStart[0], yStart[1], yStart[2]);
    finalPoint   = G4ThreeVector( yEnd[0], yEnd[1],  yEnd[2]);
    midPoint = G4ThreeVector( yMid[0], yMid[1], yMid[2]);

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

// (New) DistChord function using interpolation
G4double G4DormandPrince745::DistChord2() const
{
    // Copy the values of stages from this (original) into the Aux Stepper
    *fAuxStepper = *this;

    //Preparing for the interpolation
    fAuxStepper->SetupInterpolation(); // (fLastInitialVector, fInitialDyDx, fLastStepLength);
    //Interpolate to half step
    fAuxStepper->Interpolate( /*fLastInitialVector, fInitialDyDx, fLastStepLength,*/ 0.5, fAuxStepper->fMidVector);

    return DistLine( fLastInitialVector, fAuxStepper->fMidVector, fLastFinalVector);
}

G4double G4DormandPrince745::DistChord() const
{
    //Coefficients for halfway interpolation
    const G4double
    hf1 = 5783653.0/57600000.0 ,
    hf2 = 0. ,
    hf3 = 466123.0/1192500.0 ,
    hf4 = -41347.0/1920000.0 ,
    hf5 = 16122321.0/339200000.0 ,
    hf6 = -7117.0/20000.0,
    hf7 = 183.0/10000.0 ;

    for(int i=0; i<3; i++){
       fMidVector[i] = fLastInitialVector[i] + fLastStepLength*(
                    hf1*fInitialDyDx[i] + hf2*ak2[i] + hf3*ak3[i] + hf4*ak4[i] +
                    hf5*ak5[i] + hf6*ak6[i] + hf7*ak7[i] );
    }
     
    // Use stored values of Initial and Endpoint + new Midpoint to evaluate
    //  distance of Chord

    return DistLine( fLastInitialVector, fMidVector, fLastFinalVector);
}

//The original DistChord() function for the class
G4double  G4DormandPrince745::DistChord3() const
{
    // Do half a step using StepNoErr    
    fAuxStepper->Stepper( fLastInitialVector, fInitialDyDx, 0.5 * fLastStepLength,
                          fAuxStepper->fMidVector,  fAuxStepper->fMidError) ;
    return DistLine( fLastInitialVector, fAuxStepper->fMidVector, fLastFinalVector);
}

// The lower (4th) order interpolant given by Dormand and prince
//	"An RK 5(4) triple"
//---Ref---
//	J. R. Dormand and P. J. Prince, “Runge-Kutta triples,”
//	Computers & Mathematics with Applications, vol. 12, no. 9,
//	pp. 1007–1017, 1986.
//---------------------------

void G4DormandPrince745::SetupInterpolation_low() // const G4double *yInput, const G4double *dydx, const G4double Step)
{
    //Nothing to be done
}

void G4DormandPrince745::Interpolate_low( /* const G4double yInput[],
                                                const G4double dydx[], 
                                                const G4double Step, */
                                                G4double yOut[],
                                               G4double tau )
{
    G4double bf1, bf2, bf3, bf4, bf5, bf6, bf7;
    // Coefficients for all the seven stages.
    G4double Step = fLastStepLength;
    const G4double *dydx= fInitialDyDx;
    
    const G4int numberOfVariables= this->GetNumberOfVariables();

    // for(int i=0;i<numberOfVariables;i++) { yIn[i]=yInput[i]; }
    
    const G4double
      tau_2 = tau   * tau,
      tau_3 = tau   * tau_2,
      tau_4 = tau_2 * tau_2;
    
    bf1 = (157015080.0*tau_4 - 13107642775.0*tau_3+ 34969693132.0*tau_2- 32272833064.0*tau
           + 11282082432.0)/11282082432.0,
    bf2 = 0.0 ,
    bf3 = - 100.0*tau*(15701508.0*tau_3 - 914128567.0*tau_2 + 2074956840.0*tau
                 - 1323431896.0)/32700410799.0,
    bf4 = 25.0*tau*(94209048.0*tau_3- 1518414297.0*tau_2+ 2460397220.0*tau - 889289856.0)/5641041216.0 ,
    bf5 = -2187.0*tau*(52338360.0*tau_3 - 451824525.0*tau_2 + 687873124.0*tau - 259006536.0)/199316789632.0 ,
    bf6 =  11.0*tau*(106151040.0*tau_3- 661884105.0*tau_2 + 946554244.0*tau - 361440756.0)/2467955532.0 ,
    bf7 = tau*(1.0 - tau)*(8293050.0*tau_2 - 82437520.0*tau + 44764047.0)/ 29380423.0 ;

    //Putting together the coefficients calculated as the respective stage coefficients
    for( int i=0; i<numberOfVariables; i++){
        yOut[i] = yIn[i] + Step*tau*(bf1*dydx[i] + bf2*ak2[i] + bf3*ak3[i] + bf4*ak4[i]
                                     + bf5*ak5[i] + bf6*ak6[i] + bf7*ak7[i]  ) ;
    }
}


// Following interpolant of order 5 was given by Baker,Dormand,Gilmore, Prince :
//---Ref---
//	T. S. Baker, J. R. Dormand, J. P. Gilmore, and P. J. Prince,
//	“Continuous approximation with embedded Runge-Kutta methods,”
//	Applied Numerical Mathematics, vol. 22, no. 1, pp. 51–62, 1996.
//---------------------

// Calculating the extra stages for the interpolant :
void G4DormandPrince745::SetupInterpolation_high( /* const G4double yInput[],
                                               const G4double dydx[],
                                               const G4double Step */  ){
    
    //Coefficients for the additional stages :
    const G4double
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
    const G4double *dydx = fInitialDyDx;
    const G4double  Step = fLastStepLength;
    
    //  Saving yInput because yInput and yOut can be aliases for same array
    // for(int i=0;i<numberOfVariables;i++) { yIn[i]=yInput[i]; }
    // yTemp[7]  = yIn[7];

    //Evaluate the extra stages :
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b81*dydx[i] + b82*ak2[i] + b83*ak3[i] +
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] +
                                  b87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);				//8th Stage
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                 b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
                                 b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide(yTemp, ak9);          //9th Stage
}


// Calculating the interpolated result yOut with the coefficients
void G4DormandPrince745::Interpolate_high( /* const G4double yInput[],
                                         const G4double dydx[],
                                         const G4double Step, */
                                               G4double yOut[],
                                               G4double tau ){
    //Define the coefficients for the polynomials
    G4double bi[10][5], b[10];
    const G4int numberOfVariables = this->GetNumberOfVariables();
    const G4double *dydx = fInitialDyDx;
    // const G4double fullStep = fLastStepLength;

    // If given requestedStep in argument:
    // G4double tau = requestedStep / fLastStepLength;
    
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
    
    // for(G4int i = 0; i< numberOfVariables; i++) { yIn[i] = yInput[i]; }
    
    //    Calculating the polynomials :
#if 1    
    for(int iStage=1; iStage<=9; iStage++){
        b[iStage] = 0;
    }

    for(int j=0; j<=4; j++){
       G4double tauPower = 1.0;       
       for(int iStage=1; iStage<=9; iStage++){
            b[iStage] += bi[iStage][j]*tauPower;
       }
       tauPower *= tau;       
    }
#else    
    G4double tau0 = tau;
    
    for(int i=1; i<=9; i++){	//Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = 1.0;
        for(int j=0; j<=4; j++){
            b[i] += bi[i][j]*tau;
            tau*=tau0;
        }
    }    
#endif
    
    G4double stepLen = fLastStepLength * tau;
    for(int i=0; i<numberOfVariables; i++){		//Here i IS the cooridnate no.
        yOut[i] = yIn[i] + stepLen *(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                     b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                     b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] );
    }

}

//G4DormandPrince745::G4DormandPrince745(G4DormandPrince745& DP_Obj){
//    
//}

// Overloaded = operator
G4DormandPrince745& G4DormandPrince745::operator=(const G4DormandPrince745& right)
{
//    this->G4DormandPrince745(right.GetEquationOfMotion(),right.GetNumberOfVariables(), false);

    int noVars = right.GetNumberOfVariables();
    for(int i =0; i< noVars; i++)
    {
        this->ak2[i] = right.ak2[i];
        this->ak3[i] = right.ak3[i];
        this->ak4[i] = right.ak4[i];
        this->ak5[i] = right.ak5[i];
        this->ak6[i] = right.ak6[i];
        this->ak7[i] = right.ak7[i];
        this->ak8[i] = right.ak8[i];
        this->ak9[i] = right.ak9[i];     

        this->fInitialDyDx[i] = right.fInitialDyDx[i];
        this->fLastInitialVector[i] = right.fLastInitialVector[i];
        this->fMidVector[i] = right.fMidVector[i];
        this->fMidError[i] = right.fMidError[i];
    }
    
    this->fLastStepLength = right.fLastStepLength;
    
    return *this;
}
