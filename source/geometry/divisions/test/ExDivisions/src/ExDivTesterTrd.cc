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
// $Id: ExDivTesterTrd.cc,v 1.1 2003-11-19 18:00:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class ExDivTesterTrd Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterTrd.hh"
#include "G4Trd.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterTrd::
ExDivTesterTrd( PVType& pvtype, std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kXAxis );
  theAxis.push_back( kYAxis );
  theAxis.push_back( kZAxis );
}
 
//--------------------------------------------------------------------------
void ExDivTesterTrd::GenerateScanPoints()
{
  GenerateScanPointsAsBox();
}

//--------------------------------------------------------------------------
void ExDivTesterTrd::BuildParentSolids()
{
  theParentSolids.push_back( new G4Trd("parent_1", theWorldLengthXY, theWorldLengthXY,
                             theWorldLengthXY*0.5, theWorldLengthXY, theWorldLengthXY) );
  theParentSolids.push_back( new G4Trd("parent_2", theWorldLengthXY*0.5, theWorldLengthXY,
                             theWorldLengthXY, theWorldLengthXY, theWorldLengthXY) );
  theParentSolids.push_back( new G4Trd("parent_3", theWorldLengthXY*0.5, theWorldLengthXY,
                             theWorldLengthXY*0.5, theWorldLengthXY, theWorldLengthXY) );
}

//--------------------------------------------------------------------------
void ExDivTesterTrd::BuildChildrenSolids()
{
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );
  theWidths.push_back( 2*theWorldLengthXY / theNDiv );

  theChildSolids.push_back( new G4Trd("child_1", theWorldLengthXY*0.2,
                            theWorldLengthXY*0.2, theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY) );
  theChildSolids.push_back( new G4Trd("child_2", theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY*0.2,
                            theWorldLengthXY*0.2, theWorldLengthXY) );
  theChildSolids.push_back( new G4Trd("child_3", theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY*0.5,
                            theWorldLengthXY, theWorldLengthXY*0.2) );
}

