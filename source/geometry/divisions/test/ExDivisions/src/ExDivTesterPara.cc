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
// $Id: ExDivTesterPara.cc,v 1.1 2003-11-19 18:00:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class ExDivTesterPara Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExDivTesterPara.hh"
#include "G4Para.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4PVPlacement.hh"

//--------------------------------------------------------------------------
ExDivTesterPara::
ExDivTesterPara( PVType& pvtype, std::vector<G4String>& extraPars )
  : ExVDivTester( pvtype, extraPars )
{
  //----- Get the axis of division
  theAxis.push_back( kXAxis ); 
  theAxis.push_back( kYAxis );
  theAxis.push_back( kZAxis );
}

//--------------------------------------------------------------------------
void ExDivTesterPara::GenerateScanPoints()
{
  GenerateScanPointsAsBox();
}

//--------------------------------------------------------------------------
void ExDivTesterPara::BuildParentSolids()
{
  theParentSolids.push_back( new G4Para("parent_1", theWorldLengthXY*0.4,
                             theWorldLengthXY*0.5, theWorldLengthXY*0.6,
                             30.*deg, 45.*deg, 60.*deg ) );
  theParentSolids.push_back( new G4Para("parent_2", theWorldLengthXY*0.4,
                             theWorldLengthXY*0.5, theWorldLengthXY*0.6,
                             30.*deg, 45.*deg, 60.*deg ) );
  theParentSolids.push_back( new G4Para("parent_3", theWorldLengthXY*0.4,
                             theWorldLengthXY*0.5, theWorldLengthXY*0.6,
                             30.*deg, 45.*deg, 60.*deg ) );
}

//--------------------------------------------------------------------------
void ExDivTesterPara::BuildChildrenSolids()
{
  theChildSolids.push_back( new G4Para("child_1", theWorldLengthXY*0.4,
                            theWorldLengthXY*0.5, theWorldLengthXY*0.6/theNDiv,
                            30.*deg, 45.*deg, 60.*deg ) );
  theChildSolids.push_back( new G4Para("child_2", theWorldLengthXY*0.4,
                            theWorldLengthXY*0.5, theWorldLengthXY*0.6/theNDiv,
                            30.*deg, 45.*deg, 60.*deg ) );
  theChildSolids.push_back( new G4Para("child_3", theWorldLengthXY*0.4,
                            theWorldLengthXY*0.5, theWorldLengthXY*0.6/theNDiv,
                            30.*deg, 45.*deg, 60.*deg ) );
  theWidths.push_back( 2*((G4Para*)theParentSolids[0])->GetXHalfLength() / theNDiv );
  theWidths.push_back( 2*((G4Para*)theParentSolids[1])->GetYHalfLength() / theNDiv );
  theWidths.push_back( 2*((G4Para*)theParentSolids[2])->GetZHalfLength() / theNDiv );

}

