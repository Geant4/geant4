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
// $Id: ExDivDetectorConstruction.cc,v 1.1 2003-11-19 18:00:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "ExDivDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

#include "ExDivTesterBox.hh"
#include "ExDivTesterTubs.hh"
#include "ExDivTesterCons.hh"
#include "ExDivTesterTrd.hh"
#include "ExDivTesterPara.hh"
#include "ExDivTesterPolycone.hh"
#include "ExDivTesterPolyhedra.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExDivDetectorConstruction::
ExDivDetectorConstruction( const G4String& solidTypeStr,
                           const G4String& PVTypeStr,
                           const std::vector<G4String>& extraPars )
  : theSolidTypeStr( solidTypeStr ),
    thePVTypeStr( PVTypeStr ),
    theExtraPars( extraPars )
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExDivDetectorConstruction::~ExDivDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* ExDivDetectorConstruction::Construct()
{

  theDivTester = CreateSolidTester(theSolidTypeStr,thePVTypeStr,theExtraPars);

  //-  SolidType soltype = getSolidType( theSolidTypeStr );
  //  PVType pvtype = pvDivision;

  G4VPhysicalVolume* myTopNode = theDivTester->BuildGeometry( );
 
  theDivTester->GenerateScanPoints();

  /*  // Repeat tests but with full voxels
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4GeometryManager::GetInstance()->CloseGeometry(true);
  testG4Navigator1(myTopNode);
  testG4Navigator2(myTopNode);
  */
  
  // theDivTester->PrintParentSolid( G4cout );
  // theDivTester->PrintChildrenSolids( G4cout );
 
  // G4GeometryManager::GetInstance()->OpenGeometry();

  return myTopNode;
}

//--------------------------------------------------------------------------
ExVDivTester*
ExDivDetectorConstruction::CreateSolidTester( const G4String& stype,
                                              const G4String& thePVTypeStr,
                                              std::vector<G4String>& extraPars )
{
  PVType pvtype = getPVType( thePVTypeStr );

  ExVDivTester* theSolidTester = 0;
  if( stype == "box" ) {
   theSolidTester = new ExDivTesterBox( pvtype, extraPars );
  } else if( stype == "tubs" ) {
    theSolidTester = new ExDivTesterTubs( pvtype, extraPars );
    ExVDivTester::bDivCylindrical = 1;
  } else if( stype == "cons" ) {
    theSolidTester = new ExDivTesterCons( pvtype, extraPars );
    ExVDivTester::bDivCylindrical = 1;
  } else if( stype == "trd" ) {
    theSolidTester = new ExDivTesterTrd( pvtype, extraPars );
  } else if( stype == "para" ) {
    theSolidTester = new ExDivTesterPara( pvtype, extraPars );
  } else if( stype == "pcone" ) {
    theSolidTester = new ExDivTesterPolycone( pvtype, extraPars );
    ExVDivTester::bDivCylindrical = 1;
  } else if( stype == "phedra" ) {
    theSolidTester = new ExDivTesterPolyhedra( pvtype, extraPars );
  } else {
    G4cout << "ERROR - ExDivDetectorConstruction::CreateSolidTester()"
           << G4endl
           << "        Only the following solid types are allowed:"
           << G4endl
           << "        'tubs', 'cons', 'trd', 'para', 'pcone', 'phedra'."
           << G4endl
           << "        Wrong identifier: " << stype << G4endl;
    G4Exception("ExDivDetectorConstruction::CreateSolidTester()",
                "InvalidSetup", FatalException, "Unknown solid type.");
  }
  return theSolidTester;
}


//--------------------------------------------------------------------------
PVType ExDivDetectorConstruction::getPVType( const G4String& pvt )
{
  G4cout << pvt << G4endl;
  PVType vtype = pvPlacement;

  if( pvt == "division")
  {
    vtype = pvDivision; 
  }
  else if( pvt == "replica")
  {
    vtype = pvReplica; 
  }
  else
  {
    G4Exception("ExDivDetectorConstruction::getPVType()",
                "InvalidArgument", FatalException,
                "PV type can only be 'division' or 'replica' ");
  }
  return vtype;
}
