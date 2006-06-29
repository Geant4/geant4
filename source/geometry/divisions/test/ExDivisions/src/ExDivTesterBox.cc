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
// $Id: ExDivTesterBox.cc,v 1.3 2006-06-29 18:20:08 gunter Exp $
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
ExDivTesterBox( PVType& pvtype, PlaceType& postype,
                std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, postype, extraPars )
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

