//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExDivTesterTubs.cc,v 1.1 2003-11-19 18:00:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class ExDivTesterTubs Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterTubs.hh"
#include "G4Tubs.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterTubs::
ExDivTesterTubs( PVType& pvtype, std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kRho );
  theAxis.push_back( kPhi );
  theAxis.push_back( kZAxis );
}
 
//--------------------------------------------------------------------------
void ExDivTesterTubs::GenerateScanPoints()
{
  GenerateScanPointsAsTubs();
}

//--------------------------------------------------------------------------
void ExDivTesterTubs::BuildParentSolids()
{
  theParentSolids.push_back( new G4Tubs("parent_1", 0.5 * theWorldLengthXY,
                             theWorldLengthXY, theWorldLengthXY, theStartPhi, theDeltaPhi) );
  theParentSolids.push_back( new G4Tubs("parent_2", 0.5 * theWorldLengthXY,
                             theWorldLengthXY, theWorldLengthXY, theStartPhi, theDeltaPhi) );
  theParentSolids.push_back( new G4Tubs("parent_3", 0.5 * theWorldLengthXY,
                             theWorldLengthXY, theWorldLengthXY, theStartPhi, theDeltaPhi) );
}

//--------------------------------------------------------------------------
void ExDivTesterTubs::BuildChildrenSolids()
{
    G4Tubs* pTubs = reinterpret_cast<G4Tubs*>( theParentSolids[0] );

    theWidths.push_back( (pTubs->GetOuterRadius()-pTubs->GetInnerRadius()) / theNDiv );
    theWidths.push_back( pTubs->GetDeltaPhiAngle() / theNDiv );
    theWidths.push_back( 2*pTubs->GetZHalfLength() / theNDiv );

    theChildSolids.push_back( new G4Tubs("child_1",
                              pTubs->GetInnerRadius(),
                              pTubs->GetInnerRadius()+theWidths[0],
                              pTubs->GetZHalfLength(),
                              pTubs->GetStartPhiAngle(),
                              pTubs->GetDeltaPhiAngle() ));
    theChildSolids.push_back( new G4Tubs("child_2",
                              pTubs->GetInnerRadius(),
                              pTubs->GetOuterRadius(),
                              pTubs->GetZHalfLength(),
                              pTubs->GetStartPhiAngle(),
                              theWidths[1] ));
    theChildSolids.push_back( new G4Tubs("child_3",
                              pTubs->GetInnerRadius(),
                              pTubs->GetOuterRadius(),
                              theWidths[2]/2.,
                              pTubs->GetStartPhiAngle(),
                              pTubs->GetDeltaPhiAngle() ));
}

