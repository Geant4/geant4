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
// $Id: ExDivTesterBox.cc,v 1.1 2003-11-19 18:00:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class ExDivTesterBox Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterBox.hh"
#include "G4Box.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterBox::
ExDivTesterBox( PVType& pvtype, std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kXAxis );
  theAxis.push_back( kYAxis );
  theAxis.push_back( kZAxis );
}

 
//--------------------------------------------------------------------------
void ExDivTesterBox::GenerateScanPoints()
{
  GenerateScanPointsAsBox();
}

//--------------------------------------------------------------------------
void ExDivTesterBox::BuildParentSolids()
{
  theParentSolids.push_back( new G4Box("parent_1",
                             theWorldLengthXY, theWorldLengthXY, theWorldLengthXY) );
  theParentSolids.push_back( new G4Box("parent_2",
                             theWorldLengthXY, theWorldLengthXY, theWorldLengthXY) );
  theParentSolids.push_back( new G4Box("parent_3",
                             theWorldLengthXY, theWorldLengthXY, theWorldLengthXY) );
}

//--------------------------------------------------------------------------
void ExDivTesterBox::BuildChildrenSolids()
{
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );

  theChildSolids.push_back( new G4Box("child_1", theWorldLengthXY/theNDiv,
                            theWorldLengthXY, theWorldLengthXY) );
  theChildSolids.push_back( new G4Box("child_2", theWorldLengthXY,
                            theWorldLengthXY/theNDiv, theWorldLengthXY) );
  theChildSolids.push_back( new G4Box("child_3", theWorldLengthXY,
                            theWorldLengthXY, theWorldLengthXY/theNDiv) );
}

