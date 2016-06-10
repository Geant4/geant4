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
// $Id$

// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#include "G4PlotMessenger.hh"
#include "G4PlotParameters.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AnalysisMessengerHelper.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
//#include "G4Tokenizer.hh"

#include <iostream>
#include <sstream>
#include <vector>

using namespace G4Analysis;

//_____________________________________________________________________________
G4PlotMessenger::G4PlotMessenger(G4PlotParameters* plotParameters)
  : G4UImessenger(),
    fPlotParameters(plotParameters),
    fHelper(nullptr),
    fDirectory(nullptr),
    fSetLayoutCmd(nullptr),
    fSetDimensionsCmd(nullptr), 
    fSetStyleCmd(nullptr)
{  
  fHelper = G4Analysis::make_unique<G4AnalysisMessengerHelper>("plot");

  fDirectory = fHelper->CreateHnDirectory();

  SetStyleCmd();
  SetLayoutCmd();
  SetDimensionsCmd();
}

//_____________________________________________________________________________
G4PlotMessenger::~G4PlotMessenger()
{}

//
// private functions
//

//_____________________________________________________________________________
void G4PlotMessenger::SetStyleCmd()
{
  fSetStyleCmd = G4Analysis::make_unique<G4UIcmdWithAString>("/analysis/plot/setStyle",this);
#if defined(TOOLS_USE_FREETYPE)
  fSetStyleCmd->SetGuidance("Set plotting style from: ");
  fSetStyleCmd->SetGuidance("  ROOT_default:  ROOT style with high resolution fonts");
  fSetStyleCmd->SetGuidance("  hippodraw:     hippodraw style with high resolution fonts");
  fSetStyleCmd->SetGuidance("  inlib_default: PAW style with low resolution fonts");
  fSetStyleCmd->SetParameterName("Style", false);
#else
  fSetStyleCmd->SetGuidance("Only one plotting style is available in low resolution: ");
  fSetStyleCmd->SetGuidance("  inlib_default: PAW style with low resolution fonts");
  fSetStyleCmd->SetParameterName("Style", false);
  fSetStyleCmd->SetCandidates("inlib_default");
#endif
  fSetStyleCmd->SetCandidates(fPlotParameters->GetAvailableStyles());
  fSetStyleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
void G4PlotMessenger::SetLayoutCmd()
{
  auto columns = new G4UIparameter("columns", 'i', false);
  columns->SetGuidance("The number of columns in the page layout.");
  G4String range = "columns>=1 && columns<=";
  std::ostringstream osMaxColumns;
  osMaxColumns << fPlotParameters->GetMaxColumns();
  range.append(osMaxColumns.str());
  columns->SetParameterRange(range);
  
  auto rows = new G4UIparameter("rows", 'i', false);
  rows->SetGuidance("The number of rows in the page layout.");
  range = "rows>=1 && rows<=";
  std::ostringstream osMaxRows;
  osMaxRows << fPlotParameters->GetMaxRows();
  range.append(osMaxRows.str());
  rows->SetParameterRange(range);

  fSetLayoutCmd = G4Analysis::make_unique<G4UIcommand>("/analysis/plot/setLayout", this);
  // Guidance text:
  // Set page layout (number of columns and rows per page).
  //    Suported layouts: 
  //    columns = 1 .. maxValueAllowed
  //    rows    = 1 .. maxValueAllowed, and >= columns
  fSetLayoutCmd->SetGuidance("Set page layout (number of columns and rows per page).");
  fSetLayoutCmd->SetGuidance("   Suported layouts: ");
  G4String guidance = "  columns = 1 .. ";
  guidance.append(osMaxColumns.str());
  fSetLayoutCmd->SetGuidance(guidance);
  guidance = "  rows    = 1 .. ";
  guidance.append(osMaxRows.str());
  guidance.append(" and  >= columns");
  fSetLayoutCmd->SetGuidance(guidance);
  fSetLayoutCmd->SetParameter(columns);
  fSetLayoutCmd->SetParameter(rows);
  fSetLayoutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4PlotMessenger::SetDimensionsCmd()
{
  auto width = new G4UIparameter("width", 'i', false);
  width->SetGuidance("The page width.");
  
  auto height = new G4UIparameter("height", 'i', false);
  height->SetGuidance("The page height.");
  
  fSetDimensionsCmd = G4Analysis::make_unique<G4UIcommand>("/analysis/plot/setDimensions", this);
  fSetDimensionsCmd->SetGuidance("Set the plotter window size (width and height) in pixels.");
  fSetDimensionsCmd->SetParameter(width);
  fSetDimensionsCmd->SetParameter(height);
  fSetDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//
// public functions
//

//_____________________________________________________________________________
void G4PlotMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  // tokenize parameters in a vector
  std::vector<G4String> parameters;
  G4Analysis::Tokenize(newValues, parameters);
  // check consistency
  if ( G4int(parameters.size()) != command->GetParameterEntries() ) {
    // Should never happen but let's check anyway for consistency
    fHelper->WarnAboutParameters(command, parameters.size());
    return;
  }  

  if ( command == fSetLayoutCmd.get() ) { 
    auto counter = 0;
    auto columns = G4UIcommand::ConvertToInt(parameters[counter++]);
    auto rows = G4UIcommand::ConvertToInt(parameters[counter++]);
    fPlotParameters->SetLayout(columns, rows);
  }
  else if ( command == fSetDimensionsCmd.get() ) {
    auto counter = 0;
    auto width = G4UIcommand::ConvertToInt(parameters[counter++]);
    auto height = G4UIcommand::ConvertToInt(parameters[counter++]);
    fPlotParameters->SetDimensions(width, height);
  }
#if defined(TOOLS_USE_FREETYPE)
  else if ( command == fSetStyleCmd.get() ) {
    fPlotParameters->SetStyle(newValues);
  }
#endif
}  
