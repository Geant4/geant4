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
// $Id: ExDivTesterCons.cc,v 1.3 2006-06-29 18:20:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class ExDivTesterCons Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterCons.hh"
#include "G4Cons.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterCons::
ExDivTesterCons( PVType& pvtype, PlaceType& postype,
                 std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, postype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kRho );
  theAxis.push_back( kPhi );
  theAxis.push_back( kZAxis );
}

//--------------------------------------------------------------------------
void ExDivTesterCons::GenerateScanPoints()
{
  GenerateScanPointsAsTubs();
}

//--------------------------------------------------------------------------
void ExDivTesterCons::BuildParentSolids()
{

  theParentSolids.push_back( new G4Cons("parent_1", theWorldLengthXY*0.1,
                             theWorldLengthXY*0.5, theWorldLengthXY*0.6,
                             theWorldLengthXY, theWorldLengthXY,
                             theStartPhi, theDeltaPhi) );
  theParentSolids.push_back( new G4Cons("parent_2", theWorldLengthXY*0.1,
                             theWorldLengthXY*0.5, theWorldLengthXY*0.6,
                             theWorldLengthXY, theWorldLengthXY,
                             theStartPhi, theDeltaPhi) );
  theParentSolids.push_back( new G4Cons("parent_3", theWorldLengthXY*0.1,
                             theWorldLengthXY*0.5, theWorldLengthXY*0.6,
                             theWorldLengthXY, theWorldLengthXY,
                             theStartPhi, theDeltaPhi) );
}

//--------------------------------------------------------------------------
void ExDivTesterCons::BuildChildrenSolids()
{
  G4Cons* pCons0 = reinterpret_cast<G4Cons*>( theParentSolids[0] );  
  theWidths.push_back( (pCons0->GetOuterRadiusMinusZ()
                       -pCons0->GetInnerRadiusMinusZ()) / theNDiv );
  G4Cons* pCons1 = reinterpret_cast<G4Cons*>( theParentSolids[1] );
  theWidths.push_back( pCons1->GetDeltaPhiAngle() / theNDiv );
  G4Cons* pCons2 = reinterpret_cast<G4Cons*>( theParentSolids[2] );
  theWidths.push_back( 2*pCons2->GetZHalfLength() / theNDiv );
  
  G4double fwidthPlus = 0.;
  if( pCons0->GetInnerRadiusMinusZ() != 0. )
  {
    fwidthPlus = theWidths[0] * pCons0->GetInnerRadiusPlusZ()
                              / pCons0->GetInnerRadiusMinusZ();
  }
  theChildSolids.push_back( new G4Cons("child_1", pCons0->GetInnerRadiusMinusZ(),
                            pCons0->GetInnerRadiusMinusZ()+theWidths[0],
                            pCons0->GetInnerRadiusPlusZ(),
                            pCons0->GetInnerRadiusPlusZ()+fwidthPlus,
                            pCons0->GetZHalfLength(),
                            pCons0->GetStartPhiAngle(),
                            pCons0->GetDeltaPhiAngle() ));
  theChildSolids.push_back( new G4Cons("child_2", pCons1->GetInnerRadiusMinusZ(),
                            pCons1->GetOuterRadiusMinusZ(),
                            pCons1->GetInnerRadiusPlusZ(),
                            pCons1->GetOuterRadiusPlusZ(),
                            pCons1->GetZHalfLength(),
                            pCons1->GetStartPhiAngle(),
                            theWidths[1] ));
  theChildSolids.push_back( new G4Cons("child_3", pCons2->GetInnerRadiusMinusZ(),
                            pCons2->GetOuterRadiusMinusZ(),
                            pCons2->GetInnerRadiusPlusZ(),
                            pCons2->GetOuterRadiusPlusZ(),
                            theWidths[2]/2.,
                            pCons2->GetStartPhiAngle(),
                            pCons2->GetDeltaPhiAngle() ));
}

