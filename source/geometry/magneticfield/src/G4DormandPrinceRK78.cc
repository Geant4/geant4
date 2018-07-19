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
//  Dormand-Prince 8(7)13M non-FSAL implementation by Somnath Banerjee
//  Supported by Google as part of Google Summer of Code 2015.
//  Supervision / code review: John Apostolakis
//
// First version: 28 June 2015
//
//  Paper proposing this RK scheme:
//     Title:    "High order embedded Runge-Kutta formulae",
//     Authors:  P.J. Prince, J.R. Dormand
//     Journal of Computational and Applied Mathematics, Volume 7, Issue 1, 1981,
//       Pages 67-75, ISSN 0377-0427,
//     Reference:  DOI: 10.1016/0771-050X(81)90010-3
//       http://dx.doi.org/10.1016/0771-050X(81)90010-3.
//       (http://www.sciencedirect.com/science/article/pii/0771050X81900103)
//
//  History (condensed)
// -----------------------------
// 28 June 2015:  First version created         - S. Banerjee
//  4 July 2017:  Small fixes (Coverity issues) - J. Apostolakis

#include "G4DormandPrinceRK78.hh"
#include "G4LineSection.hh"

//Constructor
G4DormandPrinceRK78::G4DormandPrinceRK78(G4EquationOfMotion *EqRhs,
                                     G4int noIntegrationVariables,
                                     G4bool primary)
   : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
     fLastStepLength(-1.0), fAuxStepper(nullptr)     
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
    ak8 = new G4double[numberOfVariables];
    ak9 = new G4double[numberOfVariables];
    ak10 = new G4double[numberOfVariables];
    ak11 = new G4double[numberOfVariables];
    ak12 = new G4double[numberOfVariables];
    ak13 = new G4double[numberOfVariables];

    const G4int numStateVars = std::max(noIntegrationVariables, 8);
    yTemp = new G4double[numStateVars];
    yIn = new G4double[numStateVars] ;
    
    fLastInitialVector = new G4double[numStateVars] ;
    fLastFinalVector = new G4double[numStateVars] ;

    fLastDyDx = new G4double[numStateVars];
    
    fMidVector = new G4double[numStateVars];
    fMidError =  new G4double[numStateVars];

    if( primary )
    {
        fAuxStepper = new G4DormandPrinceRK78(EqRhs, numberOfVariables,
                                            !primary);
    }
}

//Destructor
G4DormandPrinceRK78::~G4DormandPrinceRK78(){
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;
    delete[] ak3;
    delete[] ak4;
    delete[] ak5;
    delete[] ak6;
    delete[] ak7;
    delete[] ak8;
    delete[] ak9;
    delete[] ak10;
    delete[] ak11;
    delete[] ak12;
    delete[] ak13;
    delete[] yTemp;
    delete[] yIn;
    
    delete[] fLastInitialVector;
    delete[] fLastFinalVector;
    delete[] fLastDyDx;
    delete[] fMidVector;
    delete[] fMidError;
    
    delete fAuxStepper;
    
}


//	The following scheme and the set of coefficients have been obtained from
//Table2. RK8(7)13M (Rational approximations
//---Ref---
//	P. J. Prince and J. R. Dormand, “High order embedded Runge-Kutta formulae,”
//	Journal of Computational and Applied Mathematics,
//	vol. 7, no. 1, pp. 67–75, Dec. 1980.
//------------------------------
//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void G4DormandPrinceRK78::Stepper(const G4double yInput[],
                                const G4double dydx[],
                                	  G4double Step,
                                      G4double yOut[],
                                	  G4double yErr[])
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    //G4double - only once
    const G4double
    
    b21 = 1.0/18,
    
    b31 = 1.0/48.0 ,
    b32 = 1.0/16.0 ,
    
    b41 = 1.0/32.0 ,
    b42 = 0.0 ,
    b43 = 3.0/32.0 ,
    
    b51 = 5.0/16.0 ,
    b52 =  0.0 ,
    b53 = -75.0/64.0 ,
    b54 =  75.0/64.0 ,
    
    b61 = 3.0/80.0 ,
    b62 = 0.0 ,
    b63 = 0.0 ,
    b64 = 3.0/16.0 ,
    b65 = 3.0/20.0 ,
    
    b71 = 29443841.0/614563906.0 ,
    b72 = 0.0 ,
    b73 = 0.0 ,
    b74 = 77736538.0/692538347.0 ,
    b75 = -28693883.0/1125000000.0 ,
    b76 = 23124283.0/1800000000.0 ,
    
    b81 = 16016141.0/946692911.0 ,
    b82 = 0.0 ,
    b83 = 0.0 ,
    b84 = 61564180.0/158732637.0 ,
    b85 = 22789713.0/633445777.0 ,
    b86 = 545815736.0/2771057229.0 ,
    b87 = -180193667.0/1043307555.0 ,
    
    b91 = 39632708.0/573591083.0 ,
    b92 = 0.0 ,
    b93 = 0.0 ,
    b94 = -433636366.0/683701615.0 ,
    b95 = -421739975.0/2616292301.0 ,
    b96 = 100302831.0/723423059.0 ,
    b97 = 790204164.0/839813087.0 ,
    b98 = 800635310.0/3783071287.0 ,
    
    b101 = 246121993.0/1340847787.0 ,
    b102 = 0.0 ,
    b103 = 0.0 ,
    b104 = -37695042795.0/15268766246.0 ,
    b105 = -309121744.0/1061227803.0 ,
    b106 =  -12992083.0/490766935.0 ,
    b107 = 6005943493.0/2108947869.0 ,
    b108 = 393006217.0/1396673457.0 ,
    b109 = 123872331.0/1001029789.0 ,
    
    b111 = -1028468189.0/846180014.0 ,
    b112 = 0.0 ,
    b113 = 0.0 ,
    b114 = 8478235783.0/508512852.0 ,
    b115 = 1311729495.0/1432422823.0 ,
    b116 = -10304129995.0/1701304382.0 ,
    b117 =  -48777925059.0/3047939560.0 ,
    b118 = 15336726248.0/1032824649.0 ,
    b119 =  -45442868181.0/3398467696.0 ,
    b1110 = 3065993473.0/597172653.0 ,
    
    b121 = 185892177.0/718116043.0 ,
    b122 = 0.0 ,
    b123 = 0.0 ,
    b124 = -3185094517.0/667107341.0 ,
    b125 = -477755414.0/1098053517.0 ,
    b126 = -703635378.0/230739211.0 ,
    b127 =  5731566787.0/1027545527.0 ,
    b128 = 5232866602.0/850066563.0 ,
    b129 = -4093664535.0/808688257.0 ,
    b1210 = 3962137247.0/1805957418.0 ,
    b1211 = 65686358.0/487910083.0 ,
    
    b131 = 403863854.0/491063109.0 ,
    b132 = 0.0 ,
    b133 = 0.0 ,
    b134 = -5068492393.0/434740067.0 ,
    b135 = -411421997.0/543043805.0 ,
    b136 = 652783627.0/914296604.0 ,
    b137 = 11173962825.0/925320556.0 ,
    b138 = -13158990841.0/6184727034.0 ,
    b139 = 3936647629.0/1978049680.0 ,
    b1310 = -160528059.0/685178525.0 ,
    b1311 = 248638103.0/1413531060.0 ,
    b1312 = 0.0 ,
    
    c1 = 14005451.0/335480064.0 ,
       // c2 = 0.0 ,
       // c3 = 0.0 ,
       // c4 = 0.0 ,
       // c5 = 0.0 ,
    c6 = -59238493.0/1068277825.0 ,
    c7 = 181606767.0/758867731.0 ,
    c8 = 561292985.0/797845732.0 ,
    c9 =  -1041891430.0/1371343529.0 ,
    c10 = 760417239.0/1151165299.0 ,
    c11 = 118820643.0/751138087.0 ,
    c12 = - 528747749.0/2220607170.0 ,
    c13 = 1.0/4.0 ,
    
    c_1 = 13451932.0/455176623.0 ,
       // c_2 = 0.0 ,
       // c_3 = 0.0 ,
       // c_4 = 0.0 ,
       // c_5 = 0.0 ,
    c_6 = -808719846.0/976000145.0 ,
    c_7 = 1757004468.0/5645159321.0 ,
    c_8 = 656045339.0/265891186.0 ,
    c_9 = -3867574721.0/1518517206.0 ,
    c_10 = 465885868.0/322736535.0 ,
    c_11 = 53011238.0/667516719.0 ,
    c_12 = 2.0/45.0 ,
    c_13 = 0.0 ,
    
    dc1 = c_1 - c1 ,
    // dc2 = c_2 - c2 ,
    // dc3 = c_3 - c3 ,
    // dc4 = c_4 - c4 ,
    // dc5 = c_5 - c5 ,
    dc6 = c_6 - c6 ,
    dc7 = c_7 - c7 ,
    dc8 = c_8 - c8 ,
    dc9 = c_9 - c9 ,
    dc10 = c_10 - c10 ,
    dc11 = c_11 - c11 ,
    dc12 = c_12 - c12 ,
    dc13 = c_13 - c13 ;
    
    //end of declaration !
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];
    //  Saving yInput because yInput and yOut can be aliases for same array
    
    for(i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    // RightHandSide(yIn, dydx) ;
    // 1st Stage - Not doing, getting passed
    
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
        yTemp[i] = yIn[i] + Step*(b71*dydx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yTemp, ak7);				//7th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b81*dydx[i] + b82*ak2[i] + b83*ak3[i] +
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] +
                                  b87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);				//8th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                  b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
                                  b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide(yTemp, ak9);          //9th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b101*dydx[i] + b102*ak2[i] + b103*ak3[i] +
                                  b104*ak4[i] + b105*ak5[i] + b106*ak6[i] +
                                  b107*ak7[i] + b108*ak8[i] + b109*ak9[i]);
    }
    RightHandSide(yTemp, ak10);          //10th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b111*dydx[i] + b112*ak2[i] + b113*ak3[i] +
                                  b114*ak4[i] + b115*ak5[i] + b116*ak6[i] +
                                  b117*ak7[i] + b118*ak8[i] + b119*ak9[i] +
                                  b1110*ak10[i]);
    }
    RightHandSide(yTemp, ak11);			//11th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b121*dydx[i] + b122*ak2[i] + b123*ak3[i] +
                                  b124*ak4[i] + b125*ak5[i] + b126*ak6[i] +
                                  b127*ak7[i] + b128*ak8[i] + b129*ak9[i] +
                                  b1210*ak10[i] + b1211*ak11[i]);
    }
    RightHandSide(yTemp, ak12);			//12th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b131*dydx[i] + b132*ak2[i] + b133*ak3[i] +
                                  b134*ak4[i] + b135*ak5[i] + b136*ak6[i] +
                                  b137*ak7[i] + b138*ak8[i] + b139*ak9[i] +
                                  b1310*ak10[i] + b1311*ak11[i] + b1312*ak12[i]);
    }
    RightHandSide(yTemp, ak13);			//13th and final Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        // Accumulate increments with proper weights
        
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + // c2 * ak2[i] + c3 * ak3[i]
                                 // + c4 * ak4[i] + c5 * ak5[i]
                                 + c6*ak6[i] +
                                 c7*ak7[i] + c8*ak8[i] +c9*ak9[i] + c10*ak10[i]
                                 + c11*ak11[i] + c12*ak12[i]  + c13*ak13[i]) ;
        
        // Estimate error as difference between 7th and
        // 8th order methods
        
        yErr[i] = Step*(dc1*dydx[i] + // dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i]
                        // + dc5*ak5[i] 
                        + dc6*ak6[i] + dc7*ak7[i] + dc8*ak8[i]
                        + dc9*ak9[i] + dc10*ak10[i] + dc11*ak11[i] + dc12*ak12[i]
                        + dc13*ak13[i] ) ;
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = dydx[i];
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}



//The DistChord() function fot the class - must define it here.
G4double  G4DormandPrinceRK78::DistChord() const
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




//------Verified------- - hackabot

