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
// class G4GeometryMessenger implementation
//
// Author: G.Cosmo, CERN
// --------------------------------------------------------------------

#include <iomanip>

#include "G4GeometryMessenger.hh"

#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4PropagatorInField.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4GeomTestVolume.hh"

//
// Constructor
//
G4GeometryMessenger::G4GeometryMessenger(G4TransportationManager* tman)
  : tmanager(tman)
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
  pchkCmd->SetGuidance( "This allows one to disable/re-enable verbosity in" );
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
  tolCmd->SetGuidance( "Define tolerance (in mm) by which overlaps reports" );
  tolCmd->SetGuidance( "should be reported. By default, all overlaps are" );
  tolCmd->SetGuidance( "reported, i.e. tolerance is set to: 0*mm." );
  tolCmd->SetParameterName( "Tolerance", true, true );
  tolCmd->SetDefaultValue( 0 );
  tolCmd->SetDefaultUnit( "mm" );
  tolCmd->SetUnitCategory( "Length" );

  verCmd = new G4UIcmdWithABool( "/geometry/test/verbosity", this );
  verCmd->SetGuidance( "Specify if running in verbosity mode or not." );
  verCmd->SetGuidance( "By default verbosity is set to ON (TRUE)." );
  verCmd->SetParameterName("verbosity",true);
  verCmd->SetDefaultValue(true);
  verCmd->AvailableForStates(G4State_Idle);

  rslCmd = new G4UIcmdWithAnInteger( "/geometry/test/resolution", this );
  rslCmd->SetGuidance( "Set the number of points on surface to be generated for" );
  rslCmd->SetGuidance( "checking overlaps." );
  rslCmd->SetParameterName("resolution",true);
  rslCmd->SetDefaultValue(10000);

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

  errCmd = new G4UIcmdWithAnInteger( "/geometry/test/maximum_errors", this );
  errCmd->SetGuidance( "Set the maximum number of overlap errors to report" );
  errCmd->SetGuidance( "for each single volume being checked." );
  errCmd->SetGuidance( "Once reached the maximum number specified, overlaps" );
  errCmd->SetGuidance( "affecting that volume further than that are simply ignored." );
  errCmd->SetParameterName("maximum_errors",true);
  errCmd->SetDefaultValue(1);

  parCmd = new G4UIcmdWithABool( "/geometry/test/check_parallel", this );
  parCmd->SetGuidance( "Check for overlaps in parallel worlds." );
  parCmd->SetGuidance( "By default, overlaps are only checked in the mass world (FALSE)." );
  parCmd->SetParameterName("check_parallel",true);
  parCmd->SetDefaultValue(true);

  recCmd = new G4UIcmdWithoutParameter( "/geometry/test/run", this );
  recCmd->SetGuidance( "Start running the recursive overlap check." );
  recCmd->SetGuidance( "Volumes are recursively asked to verify for overlaps" );
  recCmd->SetGuidance( "for points generated on the surface against their" );
  recCmd->SetGuidance( "respective mother volume and sisters at the same" );
  recCmd->SetGuidance( "level, performing for all daughters and daughters of" );
  recCmd->SetGuidance( "daughters, etc." );
  recCmd->SetGuidance( "NOTE: it may take a very long time," );
  recCmd->SetGuidance( "      depending on the geometry complexity !");
  recCmd->AvailableForStates(G4State_Idle);
}

//
// Destructor
//
G4GeometryMessenger::~G4GeometryMessenger()
{
  delete verCmd; delete recCmd; delete rslCmd;
  delete resCmd; delete rcsCmd; delete rcdCmd; delete errCmd;
  delete tolCmd;
  delete verbCmd; delete pchkCmd; delete chkCmd;
  delete geodir; delete navdir; delete testdir;
  for(auto* tvolume: tvolumes) {
      delete tvolume;
  }
}

//
// Init
//
void
G4GeometryMessenger::Init()
{
  // Create checker...
  //
  if (tvolumes.empty())
  {
    // Get all world volumes
    //
    const auto noWorlds = tmanager->GetNoWorlds();
    const auto fWorld = tmanager->GetWorldsIterator();
    for(size_t i=0;i<noWorlds;++i)
    {
        // Test the actual detector...
        //
        tvolumes.push_back(new G4GeomTestVolume(fWorld[i]));
    }
  }
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
  else if (command == pchkCmd) {
    SetPushFlag( newValues );
  }
  else if (command == tolCmd) {
    Init();
    tol = tolCmd->GetNewDoubleValue( newValues )
        * tolCmd->GetNewUnitValue( newValues );
    for(auto* tvolume: tvolumes)
    {
      tvolume->SetTolerance(tol);
    }
  }
  else if (command == verCmd) {
    Init();
    for(auto* tvolume: tvolumes)
    {
      tvolume->SetVerbosity(verCmd->GetNewBoolValue( newValues ));
    }
  }
  else if (command == rslCmd) {
    Init();
    for(auto* tvolume: tvolumes)
    {
      tvolume->SetResolution(rslCmd->GetNewIntValue( newValues ));
    }
  }
  else if (command == rcsCmd) {
    recLevel = rcsCmd->GetNewIntValue( newValues );
  }
  else if (command == rcdCmd) {
    recDepth = rcdCmd->GetNewIntValue( newValues );
  }
  else if (command == errCmd) {
    Init();
    for(auto* tvolume: tvolumes)
    {
      tvolume->SetErrorsThreshold(errCmd->GetNewIntValue( newValues ));
    }
  }
  else if (command == recCmd) {
    Init();
    G4cout << "Running geometry overlaps check..." << G4endl;
    RecursiveOverlapTest();
    G4cout << "Geometry overlaps check completed !" << G4endl;
  }
}

//
// GetCurrentValue
//
G4String
G4GeometryMessenger::GetCurrentValue( G4UIcommand* command )
{
  G4String cv = "";
  if (command == tolCmd)
  {
    cv = tolCmd->ConvertToString( tol, "mm" );
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
  if (!geomManager->IsGeometryClosed())
  {
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
  G4PropagatorInField* pField = tmanager->GetPropagatorInField();
  if (pField != nullptr)  { pField->CheckMode(mode); }
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
// Recursive Overlap Test
//
void
G4GeometryMessenger::RecursiveOverlapTest()
{
  // Close geometry if necessary
  //
  CheckGeometry();

  // Make test on single line supplied by user recursively
  //
  if (checkParallelWorlds)
  {
    for(auto* tvolume: tvolumes)
    {
      tvolume->TestRecursiveOverlap( recLevel, recDepth );
    }
  }
  else
  {
    tvolumes.front()->TestRecursiveOverlap( recLevel, recDepth );
  }
}
