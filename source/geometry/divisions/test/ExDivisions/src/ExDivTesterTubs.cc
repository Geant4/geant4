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
// $Id: ExDivTesterTubs.cc,v 1.3 2006-06-29 18:20:30 gunter Exp $
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
ExDivTesterTubs( PVType& pvtype, PlaceType& postype,
                 std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, postype, extraPars )
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

