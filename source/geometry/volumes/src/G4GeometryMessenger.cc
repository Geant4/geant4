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
// $Id: G4GeometryMessenger.cc,v 1.3 2001-11-01 18:59:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryMessenger
//
// Author: G.Cosmo, CERN.

#include "g4std/iomanip"
#include "G4GeometryMessenger.hh"

#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4GeomTestStreamLogger.hh"
#include "G4GeomTestVolume.hh"

//
// Constructor
//
G4GeometryMessenger::G4GeometryMessenger(G4TransportationManager* tman)
  : geometryOpened(true), x(0,0,0), p(0,0,1),
    extra(false), newtol(false), tol(1E-4*mm),
    tmanager(tman)
{
  geodir = new G4UIdirectory( "/geometry/" );
  geodir->SetGuidance( "Geometry control commands." );

  //
  // Geometry navigator commands
  //
  navdir = new G4UIdirectory( "/geometry/navigator/" );
  navdir->SetGuidance( "Geometry navigator control setup." );

  resCmd = new G4UIcmdWithoutParameter( "/geometry/navigator/reset", this );
  resCmd->SetGuidance( "Reset navigator and navigation history." );

  //
  // Geometry verification test commands
  //
  testdir = new G4UIdirectory( "/geometry/test/" );
  testdir->SetGuidance( "Geometry verification control setup." );
  testdir->SetGuidance( "Helps in detecting possible overlapping regions." );

  addCmd = new G4UIcmdWithABool( "/geometry/test/line_test",this );
  addCmd->SetGuidance( "Test along a single specified direction/position?" );
  addCmd->SetGuidance( "Use position and direction commands to change default." );
  addCmd->SetGuidance( "Default: position(0,0,0), direction(0,0,1)." );
  addCmd->SetParameterName( "Extra", true, true );

  tolCmd = new G4UIcmdWithADoubleAndUnit( "/geometry/test/tolerance",this );
  tolCmd->SetGuidance( "Set error tolerance value." );
  tolCmd->SetGuidance( "Default: 1E-4*mm." );
  tolCmd->SetParameterName( "Tolerance", true, true );
  tolCmd->SetDefaultUnit( "mm" );
  tolCmd->SetUnitCategory( "Surface" );

  posCmd = new G4UIcmdWith3VectorAndUnit( "/geometry/test/position", this );
  posCmd->SetGuidance( "Set starting position." );
  posCmd->SetParameterName( "X", "Y", "Z", true, true );
  posCmd->SetDefaultUnit( "cm" );

  dirCmd = new G4UIcmdWith3VectorAndUnit( "/geometry/test/direction", this );
  dirCmd->SetGuidance( "Set momentum direction." );
  dirCmd->SetGuidance( "Direction needs not to be a unit vector." );
  dirCmd->SetParameterName( "Px", "Py", "Pz", true, true );
  dirCmd->SetRange( "Px != 0 || Py != 0 || Pz != 0" );

  runCmd = new G4UIcmdWithoutParameter( "/geometry/test/run", this );
  runCmd->SetGuidance( "Start running the test." );
}


//
// Destructor
//
G4GeometryMessenger::~G4GeometryMessenger()
{
  delete addCmd;
  delete posCmd;
  delete dirCmd;
  delete runCmd;
  delete resCmd;
  delete tolCmd;
  delete geodir;
  delete navdir;
  delete testdir;
}


//
// SetNewValue
//
void
G4GeometryMessenger::SetNewValue( G4UIcommand* command, G4String newValues )
{
  if (command == resCmd) {
    ResetNavigator();
  }
  else if (command == addCmd) {
    extra = addCmd->GetNewBoolValue( newValues );
  }
  else if (command == posCmd) {
    x = posCmd->GetNew3VectorValue( newValues );
  }
  else if (command == dirCmd) {
    p = dirCmd->GetNew3VectorValue( newValues );
    if (p.mag() < DBL_MIN) {
      G4cerr << "Please specify non-zero momentum" << G4endl;
      p = G4ThreeVector(0,0,1);
    }
  }
  else if (command == tolCmd) {
    tol = tolCmd->GetNewDoubleValue( newValues );
    newtol = true;
  }
  else if (command == runCmd) {
    RunTest();
  }
}


//
// GetCurrentValue
//
G4String
G4GeometryMessenger::GetCurrentValue(G4UIcommand* command )
{
  G4String cv = "";
  if (command == addCmd) {
    cv = addCmd->ConvertToString( extra );
  }
  else if (command == posCmd) {
    cv = posCmd->ConvertToString( x, "cm" );
  }
  else if (command == tolCmd) {
    cv = tolCmd->ConvertToString( tol, "mm" );
  }
  else if (command == dirCmd) {
    cv = dirCmd->ConvertToString( p, "GeV" );
  }
  return cv;
}

//
// ResetNavigator
//
void
G4GeometryMessenger::ResetNavigator()
{
  // Close geometry if necessary
  //
  if (geometryOpened) {
    G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
    geomManager->OpenGeometry();
    geomManager->CloseGeometry(true);
    geometryOpened = false;
  }	

  // Reset navigator's state
  //
  G4ThreeVector p(0,0,0);
  tmanager->GetNavigatorForTracking()->LocateGlobalPointAndSetup(p,0,false);
}

//
// RunTest
//
void
G4GeometryMessenger::RunTest()
{
  // Create a logger to send errors/output to cout
  //
  G4GeomTestStreamLogger logger(G4std::cout);

  // Close geometry if necessary
  //
  if (geometryOpened) {
    G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
    geomManager->OpenGeometry();
    geomManager->CloseGeometry(true);
    geometryOpened = false;
  }	

  // Get the world volume
  //
  G4VPhysicalVolume* world =
    tmanager->GetNavigatorForTracking()->GetWorldVolume();

  // Test the actual detector...
  //
  G4GeomTestVolume geomTest(world, &logger);

  // Verify if error tolerance has changed
  //
  if (newtol) geomTest.SetTolerance(tol);

  // Run the test
  //
  if (extra)  // Make test on single line supplied by user
    geomTest.TestOneLine( x, p );
  else        // Apply default set of trajectories in a 3D grid pattern
    geomTest.TestCartGridXYZ();
  
  // Print out any errors, if found
  //
  geomTest.ReportErrors();

}
