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
// $Id: G4GeometryMessenger.cc,v 1.6 2010-11-10 14:06:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryMessenger
//
// Author: G.Cosmo, CERN.
//
// --------------------------------------------------------------------

#include <iomanip>
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
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4GeomTestStreamLogger.hh"
#include "G4GeomTestVolume.hh"

//
// Constructor
//
G4GeometryMessenger::G4GeometryMessenger(G4TransportationManager* tman)
  : x(0,0,0), p(0,0,1), grdRes(100,100,100), cylRes(90,50,50),
    cylfZ(0.8), cylfR(0.8), newtol(false), tol(1E-4*mm),
    recLevel(0), recDepth(-1), tmanager(tman), tlogger(0), tvolume(0)
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
  resCmd->SetGuidance( "NOTE: must be called only after kernel has been" );
  resCmd->SetGuidance( "      initialized once through the run manager!" );
  resCmd->AvailableForStates(G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger( "/geometry/navigator/verbose", this );
  verbCmd->SetGuidance( "Set run-time verbosity for the navigator." );
  verbCmd->SetGuidance(" 0 : Silent (default)");
  verbCmd->SetGuidance(" 1 : Display volume positioning and step lengths");
  verbCmd->SetGuidance(" 2 : Display step/safety info on point location");
  verbCmd->SetGuidance(" 3 : Display minimal state at -every- step");
  verbCmd->SetGuidance(" 4 : Maximum verbosity (very detailed!)");
  verbCmd->SetGuidance( "NOTE: this command has effect -only- if Geant4 has" );
  verbCmd->SetGuidance( "      been installed with the G4VERBOSE flag set!" );
  verbCmd->SetParameterName("level",true);
  verbCmd->SetDefaultValue(0);
  verbCmd->SetRange("level >=0 && level <=4");

  chkCmd = new G4UIcmdWithABool( "/geometry/navigator/check_mode", this );
  chkCmd->SetGuidance( "Set navigator in -check_mode- state." );
  chkCmd->SetGuidance( "This will cause extra checks to be applied during" );
  chkCmd->SetGuidance( "navigation. More strict and less tolerant conditions" );
  chkCmd->SetGuidance( "are applied. A run-time performance penalty may be" );
  chkCmd->SetGuidance( "observed when the -check_mode- state is activated." );
  chkCmd->SetGuidance( "NOTE: this command has effect -only- if Geant4 has" );
  chkCmd->SetGuidance( "      been installed with the G4VERBOSE flag set!" );
  chkCmd->SetParameterName("checkFlag",true);
  chkCmd->SetDefaultValue(false);
  chkCmd->AvailableForStates(G4State_Idle);

  pchkCmd = new G4UIcmdWithABool( "/geometry/navigator/push_notify", this );
  pchkCmd->SetGuidance( "Set navigator verbosity push notifications." );
  pchkCmd->SetGuidance( "This allows to disable/re-enable verbosity in" );
  pchkCmd->SetGuidance( "navigation, when tracks may get stuck and require" );
  pchkCmd->SetGuidance( "one artificial push along the direction by the" );
  pchkCmd->SetGuidance( "navigator. Notification is active by default." );
  pchkCmd->SetGuidance( "NOTE: this command has effect -only- if Geant4 has" );
  pchkCmd->SetGuidance( "      been installed with the G4VERBOSE flag set!" );
  pchkCmd->SetParameterName("pushFlag",true);
  pchkCmd->SetDefaultValue(true);
  pchkCmd->AvailableForStates(G4State_Idle);

  //
  // Geometry verification test commands
  //
  testdir = new G4UIdirectory( "/geometry/test/" );
  testdir->SetGuidance( "Geometry verification control setup." );
  testdir->SetGuidance( "Helps in detecting possible overlapping regions." );

  tolCmd = new G4UIcmdWithADoubleAndUnit( "/geometry/test/tolerance",this );
  tolCmd->SetGuidance( "Set error tolerance value." );
  tolCmd->SetGuidance( "Initial default value: 1E-4*mm." );
  tolCmd->SetParameterName( "Tolerance", true, true );
  tolCmd->SetDefaultValue( 1E-4 );
  tolCmd->SetDefaultUnit( "mm" );
  tolCmd->SetUnitCategory( "Length" );

  posCmd = new G4UIcmdWith3VectorAndUnit( "/geometry/test/position", this );
  posCmd->SetGuidance( "Set starting position for the line_test." );
  posCmd->SetParameterName( "X", "Y", "Z", true, true );
  posCmd->SetDefaultUnit( "cm" );

  dirCmd = new G4UIcmdWith3VectorAndUnit( "/geometry/test/direction", this );
  dirCmd->SetGuidance( "Set momentum direction for the line_test." );
  dirCmd->SetGuidance( "Direction needs not to be a unit vector." );
  dirCmd->SetParameterName( "Px", "Py", "Pz", true, true );
  dirCmd->SetRange( "Px != 0 || Py != 0 || Pz != 0" );

  linCmd = new G4UIcmdWithABool( "/geometry/test/line_test", this );
  linCmd->SetGuidance( "Performs test along a single specified direction/position." );
  linCmd->SetGuidance( "Use position and direction commands to change default." );
  linCmd->SetGuidance( "Initial default: position(0,0,0), direction(0,0,1)." );
  linCmd->SetGuidance( "If recursion flag is set to TRUE, the intersection checks" );
  linCmd->SetGuidance( "will be performed recursively in the geometry tree." );
  linCmd->SetParameterName("recursionFlag",true);
  linCmd->SetDefaultValue(false);
  linCmd->AvailableForStates(G4State_Idle);

  grzCmd = new G4UIcmdWith3Vector( "/geometry/test/grid_cells", this );
  grzCmd->SetGuidance( "Define resolution of grid geometry as number of cells," );
  grzCmd->SetGuidance( "specifying them for each dimension, X, Y and Z." );
  grzCmd->SetGuidance( "Will be applied to grid_test and recursive_test commands." );
  grzCmd->SetGuidance( "Initial default values: X=100, Y=100, Z=100." );
  grzCmd->SetParameterName( "X", "Y", "Z", true, true );
  grzCmd->SetDefaultValue( G4ThreeVector(100, 100, 100) );

  grdCmd = new G4UIcmdWithABool( "/geometry/test/grid_test", this );
  grdCmd->SetGuidance( "Start running the default grid test." );
  grdCmd->SetGuidance( "A grid of lines parallel to a cartesian axis is used;" );
  grdCmd->SetGuidance( "By default, only direct daughters of the mother volumes are checked." );
  grdCmd->SetGuidance( "If recursion flag is set to TRUE, the intersection checks" );
  grdCmd->SetGuidance( "will be performed recursively in the geometry tree." );
  grdCmd->SetGuidance( "NOTE: the recursion may take a very long time," );
  grdCmd->SetGuidance( "      depending on the geometry complexity !");
  grdCmd->SetParameterName("recursionFlag",true);
  grdCmd->SetDefaultValue(false);
  grdCmd->AvailableForStates(G4State_Idle);

  cyzCmd = new G4UIcmdWith3Vector( "/geometry/test/cylinder_geometry", this );
  cyzCmd->SetGuidance( "Define details of the cylinder geometry, specifying:" );
  cyzCmd->SetGuidance( "  nPhi - number of lines per Phi" );
  cyzCmd->SetGuidance( "  nZ   - number of Z points" );
  cyzCmd->SetGuidance( "  nRho - number of Rho points" );
  cyzCmd->SetGuidance( "Will be applied to the cylinder_test command." );
  cyzCmd->SetGuidance( "Initial default values: nPhi=90, nZ=50, nRho=50." );
  cyzCmd->SetParameterName( "nPhi", "nZ", "nRho", true, true );
  cyzCmd->SetDefaultValue( G4ThreeVector(90, 50, 50) );

  cfzCmd = new G4UIcmdWithADouble( "/geometry/test/cylinder_scaleZ", this );
  cfzCmd->SetGuidance( "Define the resolution of the cylinder geometry, specifying" );
  cfzCmd->SetGuidance( "the fraction scale for points along Z." );
  cfzCmd->SetGuidance( "Initial default values: fracZ=0.8" );
  cfzCmd->SetParameterName("fracZ",true);
  cfzCmd->SetDefaultValue(0.8);

  cfrCmd = new G4UIcmdWithADouble( "/geometry/test/cylinder_scaleRho", this );
  cfrCmd->SetGuidance( "Define the resolution of the cylinder geometry, specifying" );
  cfrCmd->SetGuidance( "the fraction scale for points along Rho." );
  cfrCmd->SetGuidance( "Initial default values: fracRho=0.8" );
  cfrCmd->SetParameterName("fracRho",true);
  cfrCmd->SetDefaultValue(0.8);

  cylCmd = new G4UIcmdWithABool( "/geometry/test/cylinder_test", this );
  cylCmd->SetGuidance( "Start running the cylinder test." );
  cylCmd->SetGuidance( "A set of lines in a cylindrical pattern of gradually" );
  cylCmd->SetGuidance( "increasing mesh size." );
  cylCmd->SetGuidance( "By default, only direct daughters of the mother volumes are checked." );
  cylCmd->SetGuidance( "If recursion flag is set to TRUE, the intersection checks" );
  cylCmd->SetGuidance( "will be performed recursively in the geometry tree." );
  cylCmd->SetGuidance( "NOTE: the recursion may take a very long time," );
  cylCmd->SetGuidance( "      depending on the geometry complexity !");
  cylCmd->SetParameterName("recursionFlag",true);
  cylCmd->SetDefaultValue(false);
  cylCmd->AvailableForStates(G4State_Idle);

  rcsCmd = new G4UIcmdWithAnInteger( "/geometry/test/recursion_start", this );
  rcsCmd->SetGuidance( "Set the initial level in the geometry tree for recursion." );
  rcsCmd->SetGuidance( "recursive_test will then start from the specified level." );
  rcsCmd->SetParameterName("initial_level",true);
  rcsCmd->SetDefaultValue(0);

  rcdCmd = new G4UIcmdWithAnInteger( "/geometry/test/recursion_depth", this );
  rcdCmd->SetGuidance( "Set the depth in the geometry tree for recursion." );
  rcdCmd->SetGuidance( "recursive_test will then stop after reached the specified depth." );
  rcdCmd->SetGuidance( "By default, recursion will proceed for the whole depth." );
  rcdCmd->SetParameterName("recursion_depth",true);
  rcdCmd->SetDefaultValue(-1);

  // Obsolete verification commands ...

  runCmd = new G4UIcmdWithABool( "/geometry/test/run", this );
  runCmd->SetGuidance( "Start running the default grid test." );
  runCmd->SetGuidance( "Same as the grid_test command." );
  runCmd->SetGuidance( "If recursion flag is set to TRUE, the intersection checks" );
  runCmd->SetGuidance( "will be performed recursively in the geometry tree." );
  runCmd->SetGuidance( "NOTE: the recursion may take a very long time," );
  runCmd->SetGuidance( "      depending on the geometry complexity !");
  runCmd->SetParameterName("recursionFlag",true);
  runCmd->SetDefaultValue(false);
  runCmd->AvailableForStates(G4State_Idle);

  recCmd = new G4UIcmdWithoutParameter( "/geometry/test/recursive_test", this );
  recCmd->SetGuidance( "Start running the recursive grid test." );
  recCmd->SetGuidance( "A grid of lines along a cartesian axis is recursively" );
  recCmd->SetGuidance( "to all daughters and daughters of daughters, etc." );
  recCmd->SetGuidance( "NOTE: it may take a very long time," );
  recCmd->SetGuidance( "      depending on the geometry complexity !");
  recCmd->AvailableForStates(G4State_Idle);
}

//
// Destructor
//
G4GeometryMessenger::~G4GeometryMessenger()
{
  delete linCmd; delete posCmd; delete dirCmd;
  delete grzCmd; delete grdCmd; delete recCmd; delete runCmd;
  delete rcsCmd; delete rcdCmd;
  delete cyzCmd; delete cfzCmd; delete cfrCmd; delete cylCmd;
  delete tolCmd;
  delete resCmd; delete verbCmd; delete pchkCmd; delete chkCmd;
  delete geodir; delete navdir; delete testdir;
  delete tvolume; delete tlogger;
}

//
// Init
//
void
G4GeometryMessenger::Init()
{
  // Clean old allocated loggers...
  //
  if (tlogger) delete tlogger;
  if (tvolume) delete tvolume;

  // Create a logger to send errors/output to cout
  //
  tlogger = new G4GeomTestStreamLogger(std::cout);

  // Get the world volume
  //
  G4VPhysicalVolume* world =
    tmanager->GetNavigatorForTracking()->GetWorldVolume();

  // Test the actual detector...
  //
  tvolume = new G4GeomTestVolume(world, tlogger);
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
  else if (command == verbCmd) {
    SetVerbosity( newValues );
  }
  else if (command == chkCmd) {
    SetCheckMode( newValues );
  }
  else if (command == posCmd) {
    x = posCmd->GetNew3VectorValue( newValues );
  }
  else if (command == dirCmd) {
    p = dirCmd->GetNew3VectorValue( newValues );
    if (p.mag() < DBL_MIN) {
      G4cerr << "Please specify non-zero momentum!" << G4endl;
      p = G4ThreeVector(0,0,1);
    }
  }
  else if (command == tolCmd) {
    tol = tolCmd->GetNewDoubleValue( newValues );
    newtol = true;
  }
  else if (command == linCmd) {
    Init();
    if (linCmd->GetNewBoolValue( newValues ))
      RecursiveLineTest();
    else
      LineTest();
  }
  else if ((command == grdCmd) || (command == runCmd)){
    Init();
    if (grdCmd->GetNewBoolValue( newValues ) || runCmd->GetNewBoolValue( newValues ))
      RecursiveGridTest();
    else
      GridTest();
  }
  else if (command == grzCmd) {
    grdRes = grzCmd->GetNew3VectorValue( newValues );
    if (grdRes.mag() < DBL_MIN) {
      G4cerr << "Please specify non-zero resolution!" << G4endl;
      grdRes = G4ThreeVector(100,100,100);
    }
  }
  else if (command == cyzCmd) {
    cylRes = cyzCmd->GetNew3VectorValue( newValues );
  }
  else if (command == cfzCmd) {
    cylfZ = cfzCmd->GetNewDoubleValue( newValues );
  }
  else if (command == cfrCmd) {
    cylfR = cfrCmd->GetNewDoubleValue( newValues );
  }
  else if (command == rcsCmd) {
    recLevel = rcsCmd->GetNewIntValue( newValues );
  }
  else if (command == rcdCmd) {
    recDepth = rcdCmd->GetNewIntValue( newValues );
  }
  else if (command == recCmd) {
    Init();
    RecursiveGridTest();
  }
  else if (command == cylCmd) {
    Init();
    if (cylCmd->GetNewBoolValue( newValues ))
      RecursiveCylinderTest();
    else
      CylinderTest();
  }
}

//
// GetCurrentValue
//
G4String
G4GeometryMessenger::GetCurrentValue(G4UIcommand* command )
{
  G4String cv = "";
  if (command == posCmd) {
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
// CheckGeometry
//
void
G4GeometryMessenger::CheckGeometry()
{
  // Verify that the geometry is closed
  //
  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
  if (!geomManager->IsGeometryClosed()) {
    geomManager->OpenGeometry();
    geomManager->CloseGeometry(true);
  }	
}

//
// ResetNavigator
//
void
G4GeometryMessenger::ResetNavigator()
{
  // Close geometry and reset optimisation if necessary
  //
  CheckGeometry();

  // Reset navigator's state
  //
  G4ThreeVector pt(0,0,0);
  G4Navigator* navigator = tmanager->GetNavigatorForTracking();
  navigator->LocateGlobalPointAndSetup(pt,0,false);
}

//
// Set navigator verbosity
//
void
G4GeometryMessenger::SetVerbosity(G4String input)
{
  G4int level = verbCmd->GetNewIntValue(input);
  G4Navigator* navigator = tmanager->GetNavigatorForTracking();
  navigator->SetVerboseLevel(level);
}

//
// Set navigator mode
//
void
G4GeometryMessenger::SetCheckMode(G4String input)
{
  G4bool mode = chkCmd->GetNewBoolValue(input);
  G4Navigator* navigator = tmanager->GetNavigatorForTracking();
  navigator->CheckMode(mode);
}

//
// Set navigator verbosity for push notifications
//
void
G4GeometryMessenger::SetPushFlag(G4String input)
{
  G4bool mode = pchkCmd->GetNewBoolValue(input);
  G4Navigator* navigator = tmanager->GetNavigatorForTracking();
  navigator->SetPushVerbosity(mode);
}

//
// LineTest
//
void
G4GeometryMessenger::LineTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Verify if error tolerance has changed
  //
  if (newtol) tvolume->SetTolerance(tol);

  // Make test on single line supplied by user
  //
  tvolume->TestOneLine( x, p );

  // Print out any errors, if found
  //
  tvolume->ReportErrors();
}

//
// RecursiveLineTest
//
void
G4GeometryMessenger::RecursiveLineTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Verify if error tolerance has changed
  //
  if (newtol) tvolume->SetTolerance(tol);

  // Make test on single line supplied by user recursively
  //
  tvolume->TestRecursiveLine( x, p, recLevel, recDepth );

  // Print out any errors, if found
  //
  tvolume->ReportErrors();
}

//
// GridTest
//
void
G4GeometryMessenger::GridTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Verify if error tolerance has changed
  //
  if (newtol) tvolume->SetTolerance(tol);

  // Apply set of trajectories in a 3D grid pattern
  //
  tvolume->TestCartGridXYZ( G4int(grdRes.x()),
                            G4int(grdRes.y()),
                            G4int(grdRes.z()) );
  
  // Print out any errors, if found
  //
  tvolume->ReportErrors();
}

//
// RecursiveGridTest
//
void
G4GeometryMessenger::RecursiveGridTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Verify if error tolerance has changed
  //
  if (newtol) tvolume->SetTolerance(tol);

  // Apply set of trajectories in a 3D grid pattern recursively
  //
  tvolume->TestRecursiveCartGrid( G4int(grdRes.x()),
                                  G4int(grdRes.y()),
                                  G4int(grdRes.z()),
                                  recLevel, recDepth );
  
  // Print out any errors, if found
  //
  tvolume->ReportErrors();
}

//
// CylinderTest
//
void
G4GeometryMessenger::CylinderTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Verify if error tolerance has changed
  //
  if (newtol) tvolume->SetTolerance(tol);

  // Apply default set of lines in a cylindrical pattern of gradually
  // increasing mesh size of trajectories in a 3D grid pattern
  //
  tvolume->TestCylinder(G4int(cylRes.x()),
                        G4int(cylRes.y()),
                        G4int(cylRes.z()),
                        cylfZ, cylfR, true);
  
  // Print out any errors, if found
  //
  tvolume->ReportErrors();
}

//
// RecursiveCylinderTest
//
void
G4GeometryMessenger::RecursiveCylinderTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Verify if error tolerance has changed
  //
  if (newtol) tvolume->SetTolerance(tol);

  // Apply default set of lines in a cylindrical pattern of gradually
  // increasing mesh size of trajectories in a 3D grid pattern recursively
  //
  tvolume->TestRecursiveCylinder(G4int(cylRes.x()),
                                 G4int(cylRes.y()),
                                 G4int(cylRes.z()),
                                 cylfZ, cylfR, true,
                                 recLevel, recDepth );
  
  // Print out any errors, if found
  //
  tvolume->ReportErrors();
}
