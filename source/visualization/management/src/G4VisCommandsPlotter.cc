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
// /vis/plotter commands - Guy Barrand October 2021.

#include "G4VisCommandsPlotter.hh"

#include "G4PlotterManager.hh"

#include <tools/tokenize>

#include <sstream>

#define G4warn G4cout

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/create //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterCreate::G4VisCommandPlotterCreate () {
  fpCommand = new G4UIcommand("/vis/plotter/create", this);
  fpCommand->SetGuidance("Create a named G4Plotter.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("name",'s',false);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterCreate::~G4VisCommandPlotterCreate () {delete fpCommand;}

void G4VisCommandPlotterCreate::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(newValue);
  _plotter.Reset();
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/setLayout ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterSetLayout::G4VisCommandPlotterSetLayout () {
  fpCommand = new G4UIcommand("/vis/plotter/setLayout", this);
  fpCommand->SetGuidance("Set plotter grid layout.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("columns",'i',true);
  parameter->SetDefaultValue(1);
  fpCommand->SetParameter (parameter);

  parameter = new G4UIparameter("rows",'i',true);
  parameter->SetDefaultValue(1);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterSetLayout::~G4VisCommandPlotterSetLayout () {delete fpCommand;}

void G4VisCommandPlotterSetLayout::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4String plotter;
  G4int cols,rows;
  std::istringstream is(newValue);
  is >> plotter >> cols >> rows;
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.SetLayout(cols,rows);

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/addStyle ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterAddStyle::G4VisCommandPlotterAddStyle () {
  fpCommand = new G4UIcommand("/vis/plotter/addStyle", this);
  fpCommand->SetGuidance("Add a style for a plotter.");
  fpCommand->SetGuidance("It is applied on all regions/plots of the plotter.");
  fpCommand->SetGuidance("default, ROOT_default, hippodraw are known embedded styles.");
  fpCommand->SetGuidance("reset is a keyword used to reset regions style.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("style",'s',true);
  parameter->SetDefaultValue("default");
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterAddStyle::~G4VisCommandPlotterAddStyle () {delete fpCommand;}

void G4VisCommandPlotterAddStyle::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4String plotter;
  G4String style;
  std::istringstream is(newValue);
  is >> plotter >> style;
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.AddStyle(style);
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/addRegionStyle //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterAddRegionStyle::G4VisCommandPlotterAddRegionStyle () {
  fpCommand = new G4UIcommand("/vis/plotter/addRegionStyle", this);
  fpCommand->SetGuidance("Add a style to be applied on a region.");
  fpCommand->SetGuidance("default, ROOT_default, hippodraw are known embedded styles.");
  fpCommand->SetGuidance("reset is a keyword used to reset a region style.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("region",'i',false);
//parameter->SetDefaultValue(0);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("style",'s',true);
  parameter->SetDefaultValue("default");
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterAddRegionStyle::~G4VisCommandPlotterAddRegionStyle () {delete fpCommand;}

void G4VisCommandPlotterAddRegionStyle::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String plotter;
  int region;
  G4String style;
  std::istringstream is(newValue);
  is >> plotter >> region >> style;
  if(region<0) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: bad region index " << region << "." << G4endl;
    }
    return;
  }
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.AddRegionStyle(region,style);
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/addRegionParameter //////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterAddRegionParameter::G4VisCommandPlotterAddRegionParameter () {
  fpCommand = new G4UIcommand("/vis/plotter/addRegionParameter", this);
  fpCommand->SetGuidance("Add a parameter to be set on a region.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("region",'i',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("parameter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("value",'s',false);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterAddRegionParameter::~G4VisCommandPlotterAddRegionParameter () {delete fpCommand;}

void G4VisCommandPlotterAddRegionParameter::SetNewValue (G4UIcommand* command, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  std::vector<std::string> args;
  tools::double_quotes_tokenize(newValue, args);
  if ( args.size() != command->GetParameterEntries() ) { // check consistency.
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: tokenize value problem." << G4endl;
    }
    return;
  }

  std::string plotter = args[0];
  int region = G4UIcommand::ConvertToInt(args[1].c_str());
  std::string parameter = args[2];
  std::string value = args[3];
  if(region<0) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: bad region index " << region << "." << G4endl;
    }
    return;
  }
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.AddRegionParameter(region,parameter,value);
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/clear ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterClear::G4VisCommandPlotterClear () {
  fpCommand = new G4UIcommand("/vis/plotter/clear", this);
  fpCommand->SetGuidance("Remove plottables from all regions.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterClear::~G4VisCommandPlotterClear () {delete fpCommand;}

void G4VisCommandPlotterClear::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(newValue);
  _plotter.Clear();
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/clearRegion /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterClearRegion::G4VisCommandPlotterClearRegion () {
  fpCommand = new G4UIcommand("/vis/plotter/clearRegion", this);
  fpCommand->SetGuidance("Remove plottables a region.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("region",'i',false);
//parameter->SetDefaultValue(0);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterClearRegion::~G4VisCommandPlotterClearRegion () {delete fpCommand;}

void G4VisCommandPlotterClearRegion::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String plotter;
  int region;
  std::istringstream is(newValue);
  is >> plotter >> region;
  if(region<0) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: bad region index " << region << "." << G4endl;
    }
    return;
  }
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.ClearRegion(region);
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/list ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterList::G4VisCommandPlotterList () {
  fpCommand = new G4UIcommand("/vis/plotter/list", this);
  fpCommand->SetGuidance("List plotters in the scene.");
}

G4VisCommandPlotterList::~G4VisCommandPlotterList () {delete fpCommand;}

void G4VisCommandPlotterList::SetNewValue (G4UIcommand*, G4String)
{
  G4PlotterManager::GetInstance().List();
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/add/h1 //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterAddRegionH1::G4VisCommandPlotterAddRegionH1 () {
  fpCommand = new G4UIcommand("/vis/plotter/add/h1", this);
  fpCommand->SetGuidance("Attach a 1D histogram to a plotter region.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("histo",'i',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("region",'i',true);
  parameter->SetDefaultValue(0);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterAddRegionH1::~G4VisCommandPlotterAddRegionH1 () {delete fpCommand;}

void G4VisCommandPlotterAddRegionH1::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  int hid;
  G4String plotter;
  int region;
  std::istringstream is(newValue);
  is >> hid >> plotter >> region;
  
  if(region<0) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: bad region index " << region << "." << G4endl;
    }
    return;
  }
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.AddRegionH1(region,hid);
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

/////////////////////////////////////////////////////////////////////////////
////////////// /vis/plotter/add/h2 //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
G4VisCommandPlotterAddRegionH2::G4VisCommandPlotterAddRegionH2 () {
  fpCommand = new G4UIcommand("/vis/plotter/add/h2", this);
  fpCommand->SetGuidance("Attach a 2D histogram to a plotter region.");

  G4UIparameter* parameter;
  parameter = new G4UIparameter("histo",'i',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("plotter",'s',false);
  fpCommand->SetParameter (parameter);
  
  parameter = new G4UIparameter("region",'i',true);
  parameter->SetDefaultValue(0);
  fpCommand->SetParameter (parameter);
}

G4VisCommandPlotterAddRegionH2::~G4VisCommandPlotterAddRegionH2 () {delete fpCommand;}

void G4VisCommandPlotterAddRegionH2::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  int hid;
  G4String plotter;
  int region;
  std::istringstream is(newValue);
  is >> hid >> plotter >> region;
  
  if(region<0) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<	"ERROR: bad region index " << region << "." << G4endl;
    }
    return;
  }
  
  G4Plotter& _plotter = G4PlotterManager::GetInstance().GetPlotter(plotter);
  _plotter.AddRegionH2(region,hid);
  
  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if(pScene) CheckSceneAndNotifyHandlers (pScene);
}

