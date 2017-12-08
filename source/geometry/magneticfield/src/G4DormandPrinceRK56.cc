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
//  FDormand-Prince RK 6(5) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 June 2015
//
//  G4DormandPrince745.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 26 June 2015
//
//
///////////////////////////////////////////////////////////////////////////////

#include "G4DormandPrinceRK56.hh"
#include "G4LineSection.hh"

//Constructor
G4DormandPrinceRK56::G4DormandPrinceRK56(G4EquationOfMotion *EqRhs,
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
    
    // Memory for Additional stages
    ak10 = new G4double[numberOfVariables];
    ak11 = new G4double[numberOfVariables];
    ak12 = new G4double[numberOfVariables];
    ak10_low = new G4double[numberOfVariables];
    
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
        fAuxStepper = new G4DormandPrinceRK56(EqRhs, numberOfVariables,
                                     !primary);
    }
}


//Destructor
G4DormandPrinceRK56::~G4DormandPrinceRK56(){
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
    delete[] ak10_low;    
    delete[] ak11;
    delete[] ak12;

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

void G4DormandPrinceRK56::Stepper(const G4double yInput[],
                                const G4double dydx[],
                                      G4double Step,
                                      G4double yOut[],
                                      G4double yErr[] )
                                 //   G4double nextDydx[] ) -- Output:
                                 //       endpoint DyDx ( for future FSAL version )
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    
    // Old Coefficients from
    // Table 1. RK6(5)8M
//---Ref---
//[P. J. Prince and J. R. Dormand, “High order embedded Runge-Kutta formulae,”
//    Journal of Computational and Applied Mathematics, vol. 7, no. 1, pp. 67–75,
//    Dec. 1980.
//----------------

    b21 =  1.0/10.0 ,
    
    b31 =  -2.0/81.0 ,
    b32 =  20.0/81.0 ,
    
    b41 =  615.0/1372.0 ,
    b42 =  -270.0/343.0 ,
    b43 =  1053.0/1372.0 ,
    
    b51 =  3243.0/5500.0 ,
    b52 =  -54.0/55.0 ,
    b53 = 50949.0/71500.0 ,
    b54 =  4998.0/17875.0 ,
    
    b61 = -26492.0/37125.0 ,
    b62 =  72.0/55.0 ,
    b63 =  2808.0/23375.0 ,
    b64 = -24206.0/37125.0 ,
    b65 =  338.0/459.0 ,
    
    b71 = 5561.0/2376.0 ,
    b72 =  -35.0/11.0 ,
    b73 =  -24117.0/31603.0 ,
    b74 = 899983.0/200772.0 ,
    b75 =  -5225.0/1836.0 ,
    b76 = 3925.0/4056.0 ,
    
    b81 = 465467.0/266112.0 ,
    b82 = -2945.0/1232.0 ,
    b83 = -5610201.0/14158144.0 ,
    b84 =  10513573.0/3212352.0 ,
    b85 = -424325.0/205632.0 ,
    b86 = 376225.0/454272.0 ,
    b87 = 0.0 ,
    
    c1 =  61.0/864.0 ,
    c2 =  0.0 ,
    c3 =  98415.0/321776.0 ,
    c4 =  16807.0/146016.0 ,
    c5 =  1375.0/7344.0 ,
    c6 =  1375.0/5408.0 ,
    c7 = -37.0/1120.0 ,
    c8 =  1.0/10.0 ,
    
    b91 =  61.0/864.0 ,
    b92 =  0.0 ,
    b93 =  98415.0/321776.0 ,
    b94 =  16807.0/146016.0 ,
    b95 =  1375.0/7344.0 ,
    b96 =  1375.0/5408.0 ,
    b97 = -37.0/1120.0 ,
    b98 =  1.0/10.0 ,

    dc1 =  c1  - 821.0/10800.0 ,
    dc2 =  c2 - 0.0 ,
    dc3 =  c3 - 19683.0/71825,
    dc4 =  c4 - 175273.0/912600.0 ,
    dc5 =  c5 - 395.0/3672.0 ,
    dc6 =  c6 - 785.0/2704.0 ,
    dc7 =  c7 - 3.0/50.0 ,
    dc8 =  c8 - 0.0 ,
    dc9 = 0.0;
    
    
// New Coefficients obtained from
    // Table 3 RK6(5)9FM with corrected coefficients
//---Ref---
//    T. S. Baker, J. R. Dormand, J. P. Gilmore, and P. J. Prince,
//    “Continuous approximation with embedded Runge-Kutta methods,”
//    Applied Numerical Mathematics, vol. 22, no. 1, pp. 51–62, 1996.
//------------------------------------
    
//    b21 =  1.0/9.0 ,
//    
//    b31 =  1.0/24.0 ,
//    b32 =  1.0/8.0 ,
//    
//    b41 =  1.0/16.0 ,
//    b42 =  0.0 ,
//    b43 =  3.0/16.0 ,
//    
//    b51 =  280.0/729.0 ,
//    b52 =  0.0 ,
//    b53 = -325.0/243.0 ,
//    b54 =  1100.0/729.0 ,
//    
//    b61 =  6127.0/14680.0 ,
//    b62 =  0.0 ,
//    b63 = -1077.0/734.0 ,
//    b64 =  6494.0/4037.0 ,
//    b65 = -9477.0/161480.0 ,
//    
//    b71 = -13426273320.0/14809773769.0 ,
//    b72 =  0.0 ,
//    b73 =  4192558704.0/2115681967.0 ,
//    b74 =  14334750144.0/14809773769.0 ,
//    b75 =  117092732328.0/14809773769.0 ,
//    b76 =  -361966176.0/40353607.0 ,
//    
//    b81 = -2340689.0/1901060.0 ,
//    b82 =  0.0 ,
//    b83 =  31647.0/13579.0 ,
//    b84 =  253549596.0/149518369.0 ,
//    b85 =  10559024082.0/977620105.0 ,
//    b86 = -152952.0/12173.0 ,
//    b87 = -5764801.0/186010396.0 ,
//    
//    b91 =  203.0/2880.0 ,
//    b92 =  0.0 ,
//    b93 =  0.0 ,
//    b94 =  30208.0/70785.0 ,
//    b95 =  177147.0/164560.0 ,
//    b96 = -536.0/705.0 ,
//    b97 =  1977326743.0/3619661760.0 ,
//    b98 = -259.0/720.0 ,
//    
//    
//    dc1 =  36567.0/458800.0 - b91,
//    dc2 =  0.0 - b92,
//    dc3 =  0.0 - b93,
//    dc4 =  9925984.0/27063465.0 - b94,
//    dc5 =  85382667.0/117968950.0 - b95,
//    dc6 = - 310378.0/808635.0 - b96 ,
//    dc7 =  262119736669.0/345979336560.0 - b97,
//    dc8 = - 1.0/2.0 - b98 ,
//    dc9 = -101.0/2294.0 ;
    

    //end of declaration
    
    
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
        yOut[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                  b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
                                  b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide(yOut, ak9);          //9th Stage
    
    
    for(i=0;i<numberOfVariables;i++)
    {
        
        // Estimate error as difference between 5th and
        // 6th order methods
        
        yErr[i] = Step*(  dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] 
                        + dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] + dc8*ak8[i]  
                        + dc9*ak9[i] ) ;

        // - Saving 'estimated' derivative at end-point
        // nextDydx[i] = ak9[i];
        
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = dydx[i];
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  G4DormandPrinceRK56::DistChord() const
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


// The following interpolation scheme has been obtained from
// Table 5. The RK6(5)9FM process and associated dense formula
//---Ref---
//	J. R. Dormand, M. A. Lockyer, N. E. McGorrigan, and P. J. Prince,
//	“Global error estimation with runge-kutta triples,”
//	Computers & Mathematics with Applications, vol. 18, no. 9, pp. 835–846, 1989.
//-----------------------------

// Fifth order interpolant with one extra function evaluation per step

void G4DormandPrinceRK56::SetupInterpolate_low( const G4double yInput[],
                                                      const G4double dydx[],
                                             const G4double Step ){
     const G4int numberOfVariables= this->GetNumberOfVariables();
    
    G4double
    b_101 = 33797.0/460800.0 ,
    b_102 = 0. ,
    b_103 = 0. ,
    b_104 = 27757.0/70785.0 ,
    b_105 = 7923501.0/26329600.0 ,
    b_106 = -927.0/3760.0 ,
    b_107 = -3314760575.0/23165835264.0 ,
    b_108 = 2479.0/23040.0 ,
    b_109 = 1.0/64.0 ;

    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b_101*dydx[i] + b_102*ak2[i] + b_103*ak3[i] +
                                  b_104*ak4[i] + b_105*ak5[i] + b_106*ak6[i] +
                                  b_107*ak7[i] + b_108*ak8[i] + b_109*ak9[i]);
    }
    RightHandSide(yTemp, ak10_low);          //10th Stage
}

void G4DormandPrinceRK56::Interpolate_low( const G4double yInput[],
                                    const G4double dydx[],
                                    const G4double Step,
                                    G4double yOut[],
                                    G4double tau ){
    {
        
        G4double
        bf1, bf4, bf5, bf6, bf7, bf8, bf9, bf10;
        
        G4double tau0 = tau;
        const G4int numberOfVariables= this->GetNumberOfVariables();

        for(int i=0;i<numberOfVariables;i++)
        {
            yIn[i]=yInput[i];
        }
        
        G4double
        tau_2 = tau0*tau0 ,
        tau_3 = tau0*tau_2,
        tau_4 = tau_2*tau_2;
        
        //bf2 = bf3 = 0
        bf1 = (66480.0*tau_4 - 206243.0*tau_3 + 237786.0*tau_2 - 124793.0*tau + 28800.0)/28800.0 ,
        bf4 = -16.0*tau*(45312.0*tau_3 - 125933.0*tau_2 + 119706.0*tau -40973.0)/70785.0 ,
        bf5 = -2187.0*tau*(19440.0*tau_3 - 45743.0*tau_2 + 34786.0*tau - 9293.0)/1645600.0 ,
        bf6 = tau*(12864.0*tau_3 - 30653.0*tau_2 + 23786.0*tau - 6533.0)/705.0 ,
        bf7 = -5764801.0*tau*(16464.0*tau_3 - 32797.0*tau_2 + 17574.0*tau - 1927.0)/7239323520.0 ,
        bf8 =  37.0*tau*(336.0*tau_3 - 661.0*tau_2 + 342.0*tau -31.0)/1440.0 ,
        bf9 = tau*(tau-1.0)*(16.0*tau_2 - 15.0*tau +3.0)/4.0 ,
        bf10 = 8.0*tau*(tau - 1.0)*(tau - 1.0)*(2.0*tau - 1.0) ;
        
        for( int i=0; i<numberOfVariables; i++){
            yOut[i] = yIn[i] + Step*tau*( bf1*dydx[i] +  bf4*ak4[i] + bf5*ak5[i] +
                                         bf6*ak6[i] + bf7*ak7[i] + bf8*ak8[i] +
                                         bf9*ak9[i] +  bf10*ak10_low[i] ) ;
        }
        
        
        
    }
    
}

//The following scheme and set of coefficients have been obtained from
//Table 2. Sixth order dense formula based on linear optimisation for RK6(5)9FM
//with extra stages C1O= 1/2, C11 =1/6, c12= 5/12

//---Ref---
//	T. S. Baker, J. R. Dormand, J. P. Gilmore, and P. J. Prince,
//  “Continuous approximation with embedded Runge-Kutta methods,”
//  Applied Numerical Mathematics, vol. 22, no. 1, pp. 51–62, 1996.
//--------------------


//	--- Sixth order interpolant with 3 additional stages per step ---

//Function for calculating the additional stages :
void G4DormandPrinceRK56::SetupInterpolate_high( const G4double yInput[],
                                          const G4double dydx[],
                                          const G4double Step ){
    
    //Coefficients for the additional stages :
    
    G4double
    b101 = 33797.0/460800.0 ,
    b102 = 0.0 ,
    b103 = 0.0 ,
    b104 = 27757.0/70785.0 ,
    b105 = 7923501.0/26329600.0 ,
    b106 = -927.0/3760.0 ,
    b107 = -3314760575.0/23165835264.0 ,
    b108 = 2479.0/23040.0 ,
    b109 = 1.0/64.0 ,
    
    b111 =  5843.0/76800.0 ,
    b112 =  0.0 ,
    b113 =  0.0 ,
    b114 =  464.0/2673.0 ,
    b115 =  353997.0/1196800.0 ,
    b116 = -15068.0/57105.0 ,
    b117 = -282475249.0/3644974080.0 ,
    b118 =  8678831.0/156245760.0 ,
    b119 =  116113.0/11718432.0 ,
    b1110 = -25.0/243.0 ,
    
    b121 = 15088049.0/199065600.0 ,
    b122 =  0.0 ,
    b123 =  0.0 ,
    b124 =  2.0/5.0 ,
    b125 =  92222037.0/268083200.0 ,
    b126 = -433420501.0/1528586640.0 ,
    b127 = -11549242677007.0/83630285291520.0 ,
    b128 =  2725085557.0/26167173120.0 ,
    b129 =  235429367.0/16354483200.0 ,
    b1210 = -90924917.0/1040739840.0 ,
    b1211 = -271149.0/21414400.0 ;
    

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
        yTemp[i] = yIn[i] + Step*(b101*dydx[i] + b102*ak2[i] + b103*ak3[i] +
                                  b104*ak4[i] + b105*ak5[i] + b106*ak6[i] +
                                  b107*ak7[i] + b108*ak8[i] + b109*ak9[i]);
    }
    RightHandSide(yTemp, ak10);          //10th Stage
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b111*dydx[i] + b112*ak2[i] + b113*ak3[i] +
                                  b114*ak4[i] + b115*ak5[i] + b116*ak6[i] +
                                  b117*ak7[i] + b118*ak8[i] + b119*ak9[i] +
                                  b1110*ak10[i]);
    }
    RightHandSide(yTemp, ak11);			//11th Stage
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b121*dydx[i] + b122*ak2[i] + b123*ak3[i] +
                                  b124*ak4[i] + b125*ak5[i] + b126*ak6[i] +
                                  b127*ak7[i] + b128*ak8[i] + b129*ak9[i] +
                                  b1210*ak10[i] + b1211*ak11[i]);
    }
    RightHandSide(yTemp, ak12);			//12th Stage

}



//Function to interpolate to tau(passed in) fraction of the step
void G4DormandPrinceRK56::Interpolate_high( const G4double yInput[],
                                     const G4double dydx[],
                                     const G4double Step,
                                           G4double yOut[],
                                           G4double tau )
{

    //Define the coefficients for the polynomials
    G4double bi[13][6], b[13];
    G4int numberOfVariables = this->GetNumberOfVariables();

	
    //  COEFFICIENTS OF   bi[ 1]
    bi[1][0] =  1.0 ,
    bi[1][1] = -18487.0/2880.0 ,
    bi[1][2] = 139189.0/7200.0 ,
    bi[1][3] = -53923.0/1800.0 ,
    bi[1][4] = 13811.0/600,
    bi[1][5] = -2071.0/300,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[2]
    bi[2][0] =  0.0 ,
    bi[2][1] =  0.0 ,
    bi[2][2] =  0.0 ,
    bi[2][3] =  0.0 ,
    bi[2][4] =  0.0 ,
    bi[2][5] =  0.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[3]
    bi[3][0] =  0.0 ,
    bi[3][1] =  0.0 ,
    bi[3][2] =  0.0 ,
    bi[3][3] =  0.0 ,
    bi[3][4] =  0.0 ,
    bi[3][5] =  0.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[4]
    bi[4][0] =  0.0 ,
    bi[4][1] = -30208.0/14157.0 ,
    bi[4][2] =  1147904.0/70785.0 ,
    bi[4][3] = -241664.0/5445.0 ,
    bi[4][4] =  241664.0/4719.0 ,
    bi[4][5] =  -483328.0/23595.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[5]
    bi[5][0] =  0.0 ,
    bi[5][1] = -177147.0/32912.0 ,
    bi[5][2] =  3365793.0/82280.0 ,
    bi[5][3] = -2302911.0/20570.0 ,
    bi[5][4] =  531441.0/4114.0 ,
    bi[5][5] = -531441.0/10285.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[6]
    bi[6][0] =  0.0 ,
    bi[6][1] =  536.0/141.0 ,
    bi[6][2] = -20368.0/705.0 ,
    bi[6][3] =  55744.0/705.0 ,
    bi[6][4] = -4288.0/47.0 ,
    bi[6][5] =  8576.0/235,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[7]
    bi[7][0] =  0.0 ,
    bi[7][1] = -1977326743.0/723932352.0 ,
    bi[7][2] =  37569208117.0/1809830880.0 ,
    bi[7][3] = -1977326743.0/34804440.0 ,
    bi[7][4] =  1977326743.0/30163848.0 ,
    bi[7][5] = -1977326743.0/75409620.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[8]
    bi[8][0] =  0.0 ,
    bi[8][1] =  259.0/144.0 ,
    bi[8][2] = -4921.0/360.0 ,
    bi[8][3] =  3367.0/90.0 ,
    bi[8][4] = -259.0/6.0 ,
    bi[8][5] =  259.0/15.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[9]
    bi[9][0] =  0.0 ,
    bi[9][1] =  62.0/105.0 ,
    bi[9][2] = -2381.0/525.0 ,
    bi[9][3] =  949.0/75.0 ,
    bi[9][4] = -2636.0/175.0 ,
    bi[9][5] =  1112.0/175.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[10]
    bi[10][0] =  0.0 ,
    bi[10][1] =  43.0/3.0 ,
    bi[10][2] = -1534.0/15.0 ,
    bi[10][3] =  3767.0/15.0 ,
    bi[10][4] = -1264.0/5.0 ,
    bi[10][5] =  448.0/5.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[11]
    bi[11][0] =  0.0 ,
    bi[11][1] =  63.0/5.0 ,
    bi[11][2] = -1494.0/25.0 ,
    bi[11][3] =  2907.0/25.0 ,
    bi[11][4] = -2592.0/25.0 ,
    bi[11][5] =  864.0/25.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTS OF  bi[12]
    bi[12][0] =  0.0 ,
    bi[12][1] = -576.0/35.0 ,
    bi[12][2] =  19584.0/175.0 ,
    bi[12][3] = -6336.0/25.0 ,
    bi[12][4] =  41472.0/175.0 ,
    bi[12][5] = -13824.0/175.0 ;
    //  --------------------------------------------------------

    
    for(G4int i = 0; i< numberOfVariables; i++)
        yIn[i] = yInput[i];
    
    G4double tau0 = tau;
    //    Calculating the polynomials (coefficents for the respective stages) :
    
    for(int i=1; i<=12; i++){	//Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = 1.0;
        for(int j=0; j<=5; j++){
            b[i] += bi[i][j]*tau;
            tau*=tau0;
        }
    }
    
//    Calculating the interpolation at the fraction tau of the step using the polynomial
//		coefficients and the respective stages
    
    for(int i=0; i<numberOfVariables; i++){		//Here i IS the cooridnate no.
        yOut[i] = yIn[i] + Step*tau0*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                 b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                 b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] +
                                 b[10]*ak10[i] + b[11]*ak11[i] + b[12]*ak12[i]);
    }
}

//-----Verified--------- - hackabot
