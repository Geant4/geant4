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

// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#include "G4PlotMessenger.hh"
#include "G4PlotParameters.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"

#include <vector>

using namespace G4Analysis;
using std::to_string;

//_____________________________________________________________________________
G4PlotMessenger::G4PlotMessenger(G4PlotParameters* plotParameters)
  : fPlotParameters(plotParameters)
{
  fDirectory = std::make_unique<G4UIdirectory>("/analysis/plot/");
  fDirectory->SetGuidance("Analysis batch plotting control");

  SetStyleCmd();
  SetLayoutCmd();
  SetDimensionsCmd();
}

//_____________________________________________________________________________
G4PlotMessenger::~G4PlotMessenger() = default;

//
// private functions
//

//_____________________________________________________________________________
void G4PlotMessenger::AddIntParameter(
  G4UIcommand& command, G4String name, G4String guidance, G4String range)
{
  auto param = new G4UIparameter(name.c_str(), 'i', false);
  param->SetGuidance(guidance.c_str());
  if (! range.empty()) {
    param->SetParameterRange(range);
  }
  command.SetParameter(param);
}

//_____________________________________________________________________________
void G4PlotMessenger::SetStyleCmd()
{
  G4String guidance;
  G4String candidates;
#if defined(TOOLS_USE_FREETYPE)
  guidance =
    "Set plotting style from: \n"
    "  ROOT_default:  ROOT style with high resolution fonts\n"
    "  hippodraw:     hippodraw style with high resolution fonts\n"
    "  inlib_default: PAW style with low resolution fonts";
  candidates =
    "ROOT_default hippodraw inlib_default";
#else
  guidance =
    "Only one plotting style is available in low resolution: \n"
    "  inlib_default: PAW style with low resolution fonts";
  candidates =
    "inlib_default";
#endif

  fSetStyleCmd = CreateCommand<G4UIcmdWithAString>("setStyle", guidance);
  fSetStyleCmd->SetParameterName("Style", false);
  fSetStyleCmd->SetCandidates("inlib_default");
}

//_____________________________________________________________________________
void G4PlotMessenger::SetLayoutCmd()
{
  fSetLayoutCmd = CreateCommand<G4UIcommand>(
    "setLayout",
    "Set page layout (number of columns and rows per page).\n"
    "   Supported layouts:\n"
    "   columns = 1 .. maxValueAllowed\n"
    "   rows    = 1 .. maxValueAllowed, and >= columns\"");

  AddIntParameter(*fSetLayoutCmd, "columns",
     "The number of columns in the page layout.",
     "columns>=1 && columns<=" + std::to_string(fPlotParameters->GetMaxColumns()));
  AddIntParameter(*fSetLayoutCmd, "rows",
    "The number of rows in the page layout.",
    "rows>=1 && rows<=" + std::to_string(fPlotParameters->GetMaxRows()));
}

//_____________________________________________________________________________
void G4PlotMessenger::SetDimensionsCmd()
{
  fSetDimensionsCmd = CreateCommand<G4UIcommand>(
    "setDimensions",
    "Set the plotter window size (width and height) in pixels.");

  AddIntParameter(*fSetDimensionsCmd, "width", "The page width.");
  AddIntParameter(*fSetDimensionsCmd, "height", "The page height.");
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
  if ( parameters.size() != command->GetParameterEntries() ) {
    // Should never happen but let's check anyway for consistency
    G4Analysis::Warn(
      "Got wrong number of \"" + command->GetCommandName() +
      "\" parameters: " + to_string(parameters.size()) +
      " instead of " + to_string(command->GetParameterEntries()) + " expected",
      fkClass, "WarnAboutParameters");
    return;
  }

  auto counter = 0;
  if ( command == fSetLayoutCmd.get() ) {
    auto columns = G4UIcommand::ConvertToInt(parameters[counter++]);
    auto rows = G4UIcommand::ConvertToInt(parameters[counter++]);
    fPlotParameters->SetLayout(columns, rows);
    return;
  }

  if ( command == fSetDimensionsCmd.get() ) {
    auto width = G4UIcommand::ConvertToInt(parameters[counter++]);
    auto height = G4UIcommand::ConvertToInt(parameters[counter++]);
    fPlotParameters->SetDimensions(width, height);
    return;
  }

  if ( command == fSetStyleCmd.get() ) {
    fPlotParameters->SetStyle(newValues);
    return;
  }
}
