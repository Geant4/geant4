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
//  Bogacki-Shampine - 8 - 5(4) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 May 2015
//
//  History
// -----------------------------
//  Created by Somnath on 26 May 2015
// 
///////////////////////////////////////////////////////////////////////////////
//  Renamed to G4 standard naming
//  Plan is that this source file / class will be merged with the updated
//  BogackiShampine45 class, which contains improvements (May 2016) 
//               J. Apostolakis,  31 May 2016
///////////////////////////////////////////////////////////////////////////////
//
//
//This is the source file of BogackiShampine45 class containing the
//definition of the stepper() method that evaluates one step in
//field propagation.
//The Butcher table of the Bogacki-Shampine-8-4-5 method is as follows :
//
//0   |
//1/6 | 1/6
//2/9 | 2/27          4/27
//3/7 | 183/1372  -162/343   1053/1372
//2/3 | 68/297      -4/11      42/143 1960/3861
//3/4 | 597/22528     81/352     63099/585728     58653/366080   4617/20480
//1   | 174197/959244 -30942/79937 8152137/19744439 666106/1039181 -29421/29068 482048/414219
//1   | 587/8064           0       4440339/15491840 24353/124800   387/44800     2152/5985   7267/94080
//-------------------------------------------------------------------------------------------------------------------
//      587/8064           0       4440339/15491840 24353/124800    387/44800     2152/5985   7267/94080       0
//      2479/34992         0          123/416       612941/3411720  43/1440       2272/6561  79937/1113912  3293/556956

#include <cassert>

#include "G4FSALBogackiShampine45.hh"
#include "G4LineSection.hh"

G4bool   G4FSALBogackiShampine45::fPreparedConstants= false;
G4double G4FSALBogackiShampine45::bi[12][7];

//Constructor
G4FSALBogackiShampine45::G4FSALBogackiShampine45(G4EquationOfMotion *EqRhs,
                                     G4int noIntegrationVariables,
                                     G4bool primary)
   : G4VFSALIntegrationStepper(EqRhs, noIntegrationVariables),
     fLastStepLength( -1.0 ), fAuxStepper( nullptr )
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
    DyDx = new G4double[numberOfVariables];
    
    assert ( GetNumberOfStateVariables() >= 8 );
    const G4int numStateVars = std::max(noIntegrationVariables,
                                        GetNumberOfStateVariables() );  

    // Must ensure space extra 'state' variables exists - i.e. yIn[7]
    yTemp = new G4double[numStateVars];
    yIn = new G4double[numStateVars] ;
    
    fLastInitialVector = new G4double[numStateVars] ;
    fLastFinalVector = new G4double[numStateVars] ;
    fLastDyDx = new G4double[numberOfVariables];  // Only derivatives
    
    fMidVector = new G4double[numStateVars];
    fMidError =  new G4double[numStateVars];
    
    pseudoDydx_for_DistChord = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    if( primary )
    {
        fAuxStepper = new G4FSALBogackiShampine45(EqRhs, numberOfVariables,
                                            !primary);
    }
    if( ! fPreparedConstants )
       PrepareConstants();
}


//Destructor
G4FSALBogackiShampine45::~G4FSALBogackiShampine45(){
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
    delete[] DyDx;    
    delete[] yTemp;
    delete[] yIn;
    
    delete[] fLastInitialVector;
    delete[] fLastFinalVector;
    delete[] fLastDyDx;
    delete[] fMidVector;
    delete[] fMidError;
    
    delete fAuxStepper;
    
    delete[] pseudoDydx_for_DistChord;
}


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void G4FSALBogackiShampine45::Stepper(const G4double yInput[],
                                const G4double dydx[],
                                G4double Step,
                                G4double yOut[],
                                G4double yErr[],
                                G4double nextDydx[])
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    
    b21 = 1.0/6.0 ,
    b31 = 2.0/27.0 , b32 = 4.0/27.0,
    
    b41 = 183.0/1372.0 , b42 = -162.0/343.0, b43 = 1053.0/1372.0,
    
    b51 = 68.0/297.0, b52 = -4.0/11.0,
    b53 = 42.0/143.0, b54 = 1960.0/3861.0,
    
    b61 = 597.0/22528.0, b62 = 81.0/352.0,
    b63 = 63099.0/585728.0, b64 = 58653.0/366080.0,
    b65 = 4617.0/20480.0,
    
    b71 = 174197.0/959244.0, b72 = -30942.0/79937.0,
    b73 = 8152137.0/19744439.0, b74 = 666106.0/1039181.0,
    b75 = -29421.0/29068.0,  b76 = 482048.0/414219.0,
    
    b81 = 587.0/8064.0,  b82 = 0.0,
    b83 = 4440339.0/15491840.0,  b84 = 24353.0/124800.0,
    b85 = 387.0/44800.0, b86 = 2152.0/5985.0,
    b87 = 7267.0/94080.0,
    
    
//    c1 = 2479.0/34992.0,
//    c2 = 0.0,
//    c3 = 123.0/416.0,
//    c4 = 612941.0/3411720.0,
//    c5 = 43.0/1440.0,
//    c6 = 2272.0/6561.0,
//    c7 = 79937.0/1113912.0,
//    c8 = 3293.0/556956.0,
    
    //For the embedded higher order method only the difference of values
    // taken and is used directly later instead of defining the last row
    // of butcher table in a separate set of variables and taking the
    // difference there

    
    dc1 = b81 - 2479.0/34992.0 ,
    dc2 = 0.0,
    dc3 = b83 - 123.0/416.0 ,
    dc4 = b84 - 612941.0/3411720.0,
    dc5 = b85 - 43.0/1440.0,
    dc6 = b86 - 2272.0/6561.0,
    dc7 = b87 - 79937.0/1113912.0,
    dc8 = -3293.0/556956.0;   //end of declaration
    
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];
    //  Saving yInput because yInput and yOut can be aliases for same array
    
    for(i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
        DyDx[i] = dydx[i];
    }
    
    
    // RightHandSide(yIn, dydx) ;
    // 1st Step - Not doing, getting passed
    
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
        yTemp[i] = yIn[i] + Step*(b41*DyDx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    RightHandSide(yTemp, ak4) ;              // 4th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b51*DyDx[i] + b52*ak2[i] + b53*ak3[i] +
                                  b54*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5) ;              // 5th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*DyDx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b71*DyDx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yTemp, ak7);              //7th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b81*DyDx[i] + b82*ak2[i] + b83*ak3[i] +
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] +
                                  b87*ak7[i]);
    }
    RightHandSide(yOut, ak8);               //8th Step - Final one Using FSAL

    
    for(i=0;i<numberOfVariables;i++)
    {
        
        yErr[i] = Step*(dc1*DyDx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] + dc8*ak8[i]) ;
        
        
        //FSAL stepper : Must pass the last DyDx for the next step, here ak8
        nextDydx[i] = ak8[i];
        
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = DyDx[i];
        
    }
    
    fLastStepLength = Step;
    
    return ;
}
//
//G4double* G4FSALBogackiShampine45::getLastDydx(){
//    return ak8;
//}

//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  G4FSALBogackiShampine45::DistChord() const
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

// ---------------------------------------------------------------------------------------

void G4FSALBogackiShampine45::PrepareConstants()
{
    //  --------------------------------------------------------
    //  COEFFICIENTS FOR INTERPOLANT  bi  WITH  11  STAGES
    //  --------------------------------------------------------
    
    // Initialise all values of G4double bi[12][7] 
    for(int i=1; i<12; i++){
       for(int j=1; j<7; j++){
          bi[i][j] = 0.0 ;
       }
    }
  
    bi[1][6] = -12134338393.0/1050809760.0 ,
    bi[1][5] = -1620741229.0/50038560.0 ,
    bi[1][4] = -2048058893.0/59875200.0 ,
    bi[1][3] = -87098480009.0/5254048800.0 ,
    bi[1][2] = -11513270273.0/3502699200.0 ,
    //
    bi[3][6] = -33197340367.0/1218433216.0 ,
    bi[3][5] = -539868024987.0/6092166080.0 ,
    bi[3][4] = -39991188681.0/374902528.0 ,
    bi[3][3] = -69509738227.0/1218433216.0 ,
    bi[3][2] = -29327744613.0/2436866432.0 ,
    //
    bi[4][6] = -284800997201.0/19905339168.0 ,
    bi[4][5] = -7896875450471.0/165877826400.0 ,
    bi[4][4] = -333945812879.0/5671036800.0 ,
    bi[4][3] = -16209923456237.0/497633479200.0 ,
    bi[4][2] = -2382590741699.0/331755652800.0 ,
    //
    bi[5][6] = -540919.0/741312.0 ,
    bi[5][5] = -103626067.0/43243200.0 ,
    bi[5][4] = -633779.0/211200.0 ,
    bi[5][3] = -32406787.0/18532800.0 ,
    bi[5][2] = -36591193.0/86486400.0 ,
    //
    bi[6][6] = 7157998304.0/374350977.0 ,
    bi[6][5] = 30405842464.0/623918295.0 ,
    bi[6][4] = 183022264.0/5332635.0 ,
    bi[6][3] = -3357024032.0/1871754885.0 ,
    bi[6][2] = -611586736.0/89131185.0 ,
    //
    bi[7][6] = -138073.0/9408.0 ,
    bi[7][5] = -719433.0/15680.0 ,
    bi[7][4] = -1620541.0/31360.0 ,
    bi[7][3] = -385151.0/15680.0 ,
    bi[7][2] = -65403.0/15680.0 ,
    //
    bi[8][6] = 1245.0/64.0 ,
    bi[8][5] = 3991.0/64.0 ,
    bi[8][4] = 4715.0/64.0 ,
    bi[8][3] = 2501.0/64.0 ,
    bi[8][2] = 149.0/16.0 ,
    bi[8][1] = 1.0 ,
    //
    bi[9][6] = 55.0/3.0 ,
    bi[9][5] = 71.0 ,
    bi[9][4] = 103.0 ,
    bi[9][3] = 199.0/3.0 ,
    bi[9][2] = 16.0 ,
    //
    bi[10][6] = -1774004627.0/75810735.0 ,
    bi[10][5] = -1774004627.0/25270245.0 ,
    bi[10][4] = -26477681.0/359975.0 ,
    bi[10][3] = -11411880511.0/379053675.0 ,
    bi[10][2] = -423642896.0/126351225.0 ,
    //
    bi[11][6] = 35.0 ,
    bi[11][5] = 105.0 ,
    bi[11][4] = 117.0 ,
    bi[11][3] = 59.0 ,
    bi[11][2] = 12.0 ;
}

// ---------------------------------------------------------------------------------------

void G4FSALBogackiShampine45::interpolate( const G4double yInput[],
                             const G4double dydx[],
                             G4double yOut[],
                             G4double Step,
                             G4double tau
                             )
{
    const G4double
    a91 = 455.0/6144.0 ,
    a92 = 0.0 ,
    a93 = 10256301.0/35409920.0 ,
    a94 = 2307361.0/17971200.0 ,
    a95 = -387.0/102400.0 ,
    a96 = 73.0/5130.0 ,
    a97 = -7267.0/215040.0 ,
    a98 = 1.0/32.0 ,

    a101 = -837888343715.0/13176988637184.0 ,
    a102 = 30409415.0/52955362.0 ,
    a103 = -48321525963.0/759168069632.0 ,
    a104 = 8530738453321.0/197654829557760.0 ,
    a105 = 1361640523001.0/1626788720640.0 ,
    a106 = -13143060689.0/38604458898.0 ,
    a107 = 18700221969.0/379584034816.0 ,
    a108 = -5831595.0/847285792.0 ,
    a109 = -5183640.0/26477681.0 ,

    a111 = 98719073263.0/1551965184000.0 ,
    a112 = 1307.0/123552.0 ,
    a113 = 4632066559387.0/70181753241600.0 ,
    a114 = 7828594302389.0/382182512025600.0 ,
    a115 = 40763687.0/11070259200.0 ,
    a116 = 34872732407.0/224610586200.0 ,
    a117 = -2561897.0/30105600.0 ,
    a118 =  1.0/10.0 ,
    a119 = -1.0/10.0 ,
    a1110 = -1403317093.0/11371610250.0 ;

    const G4int numberOfVariables= this->GetNumberOfVariables();

    //  Saving yInput because yInput and yOut can be aliases for same array
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];

    //    calculating extra stages
    for(int i=0; i<numberOfVariables; i++){
        yTemp[i] = yIn[i] + Step*(a91*dydx[i] + a92*ak2[i] + a93*ak3[i] +
                                  a94*ak4[i] + a95*ak5[i] + a96*ak6[i] +
                                  a97*ak7[i] + a98*ak8[i] );
    }
    
    RightHandSide(yTemp, ak9);
    
    for(int i=0; i<numberOfVariables; i++){
        yTemp[i] = yIn[i] + Step*(a101*dydx[i] + a102*ak2[i] + a103*ak3[i] +
                                  a104*ak4[i] + a105*ak5[i] + a106*ak6[i] +
                                  a107*ak7[i] + a108*ak8[i] + a109*ak9[i] );
    }
    
    RightHandSide(yTemp, ak10);

    for(int i=0; i<numberOfVariables; i++){
        yTemp[i] = yIn[i] + Step*(a111*dydx[i] + a112*ak2[i] + a113*ak3[i] +
                                  a114*ak4[i] + a115*ak5[i] + a116*ak6[i] +
                                  a117*ak7[i] + a118*ak8[i] + a119*ak9[i] +
                                  a1110*ak10[i] );
    }
    
    RightHandSide(yTemp, ak11);
    
    G4double tau0 = tau;
    //    Calculating the polynomials :
    for(int i=1; i<=11; i++){   //Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0.0;
        tau = tau0;
        for(int j=1; j<=6; j++){
            b[i] += bi[i][j]*tau;
            tau*=tau0;
        }
    }
    
    for(int i=0; i<numberOfVariables; i++){
        yOut[i] = yIn[i] + Step*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                 b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                 b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] +
                                 b[10]*ak10[i] + b[11]*ak11[i] );
    }  
}




