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
// G4DriverReporter
//
// Implementation
//
//
// Authors: J.Apostolakis        -   January/March 2020
// -------------------------------------------------------------------

#include "G4DriverReporter.hh"

// ---------------------------------------------------------------------------

void G4DriverReporter::PrintStatus( const G4double* StartArr,  
                                          G4double  xstart,
                                    const G4double* CurrentArr, 
                                          G4double  xcurrent,
                                          G4double  requestStep, 
                                    unsigned int    subStepNo,
                                    unsigned int    noIntegrationVariables
   )
  // Potentially add as arguments:  
  //                                 <dydx>           - as Initial Force
  //                                 stepTaken(hdid)  - last step taken
  //                                 nextStep (hnext) - proposal for size
{
   G4FieldTrack  StartFT(G4ThreeVector(0,0,0),
                 G4ThreeVector(0,0,0), 0., 0., 0., 0. );
   G4FieldTrack  CurrentFT (StartFT);

   StartFT.LoadFromArray( StartArr, noIntegrationVariables); 
   StartFT.SetCurveLength( xstart);
   CurrentFT.LoadFromArray( CurrentArr, noIntegrationVariables); 
   CurrentFT.SetCurveLength( xcurrent );

   PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
}

// ---------------------------------------------------------------------------
const G4int noPrecision = 8;
const G4int prec7= noPrecision+2;
const G4int prec8= noPrecision+3;    
const G4int prec9= noPrecision+4;

void G4DriverReporter::PrintStatus(const G4FieldTrack& StartFT,
                                   const G4FieldTrack& CurrentFT, 
                                         G4double      requestStep, 
                                   unsigned int        subStepNo)
{
    G4int verboseLevel= 2; // fVerboseLevel;
    G4long oldPrec= G4cout.precision(noPrecision);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);
        
    const G4ThreeVector StartPosition=       StartFT.GetPosition();
    const G4ThreeVector StartUnitVelocity=   StartFT.GetMomentumDir();
    const G4ThreeVector CurrentPosition=     CurrentFT.GetPosition();
    const G4ThreeVector CurrentUnitVelocity= CurrentFT.GetMomentumDir();

    G4double  DotStartCurrentVeloc= StartUnitVelocity.dot(CurrentUnitVelocity);

    G4double step_len= CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
    G4double subStepSize = step_len;
     
    if( (subStepNo <= 1) || (verboseLevel > 3) )
    {
       subStepNo = - subStepNo;        // To allow printing banner

       G4cout << "------------------------------------------------------------------"
              << G4endl;
       G4cout << std::setw( 6)  << " " << std::setw( 25)
              << " G4DriverReporter: Current Position  and  Direction" << " "
              << G4endl; 
       G4cout << std::setw( 5) << "Step#" << " "
              << std::setw( prec7) << "s-curve" << " "
              << std::setw( prec9) << "X(mm)" << " "
              << std::setw( prec9) << "Y(mm)" << " "  
              << std::setw( prec9) << "Z(mm)" << " "
              << std::setw( prec8) << " N_x " << " "
              << std::setw( prec8) << " N_y " << " "
              << std::setw( prec8) << " N_z " << " "
              << std::setw( 6) << " N^2-1 " << " "
              << std::setw(10) << " N(0).N " << " "
              << std::setw( 7) << "KinEner " << " "
              << std::setw(12) << "Track-l" << " "   // Add the Sub-step ??
              << std::setw(12) << "Step-len" << " " 
              << std::setw(12) << "Step-len" << " " 
              << std::setw( 9) << "ReqStep" << " "  
              << G4endl;
    }

    G4cout.precision(noPrecision);
    
    if( (subStepNo <= 0) )
    {
      PrintStat_Aux( StartFT,  requestStep, 0., 
                       0,        0.0,         1.0);
    }

    // if( verboseLevel <= 3 )
    {
      G4cout.precision(noPrecision);
      PrintStat_Aux( CurrentFT, requestStep, step_len, 
                     subStepNo, subStepSize, DotStartCurrentVeloc );
    }
    G4cout << "------------------------------------------------------------------"
           << G4endl;
    G4cout.precision(oldPrec);
}

// ---------------------------------------------------------------------------

void G4DriverReporter::PrintStat_Aux(const G4FieldTrack& aFieldTrack,
                                          G4double      requestStep, 
                                          G4double      step_len,
                                          G4int         subStepNo,
                                          G4double      subStepSize,
                                          G4double      dotVeloc_StartCurr)
{
    const G4ThreeVector Position = aFieldTrack.GetPosition();
    const G4ThreeVector UnitVelocity = aFieldTrack.GetMomentumDir();

    G4long oldprec= G4cout.precision(noPrecision);
    
    if( subStepNo >= 0)
    {
       G4cout << std::setw( 5) << subStepNo << " ";
    }
    else
    {
       G4cout << std::setw( 5) << "Start" << " ";
    }
    G4double curveLen= aFieldTrack.GetCurveLength();
    G4cout << std::setw( 7) << curveLen;
    // G4cout.precision(noPrecision);
    G4cout << std::setw( prec9) << Position.x() << " "
           << std::setw( prec9) << Position.y() << " "
           << std::setw( prec9) << Position.z() << " "
           << std::setw( prec8) << UnitVelocity.x() << " "
           << std::setw( prec8) << UnitVelocity.y() << " "
           << std::setw( prec8) << UnitVelocity.z() << " ";
    G4cout.precision(3);
    G4double unitMagDif = UnitVelocity.mag2()-1.0;
    if( std::fabs( unitMagDif ) < 1.0e-15 ) { unitMagDif = 0.0; }        
    G4cout << std::setw( 8) << unitMagDif << " ";
    G4cout.precision(6);
    G4cout << std::setw(10) << dotVeloc_StartCurr << " ";
    G4cout.precision(oldprec);
    G4cout << std::setw( prec7) << aFieldTrack.GetKineticEnergy();
    G4cout << std::setw(12) << step_len << " ";

    static G4ThreadLocal G4double oldCurveLength = 0.0;
    static G4ThreadLocal G4double oldSubStepLength = 0.0;
    static G4ThreadLocal G4int oldSubStepNo = -1;

    G4double subStep_len = 0.0;
    if( curveLen > oldCurveLength )
    {
      subStep_len= curveLen - oldCurveLength;
    }
    else if (subStepNo == oldSubStepNo)
    {
      subStep_len= oldSubStepLength;
    }
    oldCurveLength= curveLen;
    oldSubStepLength= subStep_len;

    G4cout << std::setw(12) << subStep_len << " "; 
    G4cout << std::setw(12) << subStepSize << " "; 
    if( requestStep != -1.0 )
    {
      G4cout << std::setw( prec9) << requestStep << " ";
    }
    else
    {
       G4cout << std::setw( prec9) << " InitialStep " << " ";
    }
    G4cout << G4endl;
}
