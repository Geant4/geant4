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

// Author: Ivana Hrivnacova, 05/05/2015  (ivana@ipno.in2p3.fr)
//
// This messenger class is a generalization of the HistoMessenger class,
// originally developed for the extended/electromagnetic examples
// by Michel Maire (michel.maire@lapp.in2p3.fr)

#include "G4AnalysisMessengerHelper.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace G4Analysis;

namespace {

//_____________________________________________________________________________
G4String ObjectType(const G4String& hnType)
{
  G4String first = hnType.substr(0,1);
  if (first == "h") {
    return "Histogram";
  } else if (first == "p") {
    return "Profile";    
  } else {
    // other possibilitied not handled
    return "";
  }
}

//_____________________________________________________________________________
void Replace(std::string& str, const std::string& from, const std::string& to) {
  // Replace all occurences of from string
  if (from.empty()) return;
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }  
} 

}                

//_____________________________________________________________________________
G4AnalysisMessengerHelper::G4AnalysisMessengerHelper(const G4String& hnType)
  : fHnType(hnType)
{}

//_____________________________________________________________________________
G4AnalysisMessengerHelper::~G4AnalysisMessengerHelper()
{}

//
// private functions
//

//_____________________________________________________________________________
G4String G4AnalysisMessengerHelper::Update(const G4String& str, const G4String& axis) const
{
  G4String newStr(str);

  // Hn, Pn
  G4String upperHnType(str);
  upperHnType.toUpper();
  Replace(newStr, "UHNTYPE_", upperHnType);

  // hn, pn
  Replace(newStr, "HNTYPE_", fHnType);
  
  // n = 1,2,3
  G4String second = fHnType.substr(1,1);
  Replace(newStr, "NDIM_", second);
  
  // histogram, profile
  G4String lowerObjectType(ObjectType(fHnType));
  lowerObjectType.toLower();
  Replace(newStr, "LOBJECT", lowerObjectType);
  
  // Histogram, Profile
  Replace(newStr, "OBJECT", ObjectType(fHnType));
  
  // X, Y, Z
  G4String upperAxis(axis);
  upperAxis.toUpper();
  Replace(newStr, "UAXIS", upperAxis);
  
  // x, y, z
  Replace(newStr, "AXIS", axis);

  // return result
  return newStr;
}

//
// public functions
//

//_____________________________________________________________________________
std::unique_ptr<G4UIdirectory> 
G4AnalysisMessengerHelper::CreateHnDirectory() const
{
  std::unique_ptr<G4UIdirectory> directory(new G4UIdirectory(Update("/analysis/HNTYPE_/")));
  directory->SetGuidance(Update("NDIM_D LOBJECT control"));
  return directory;
}

//_____________________________________________________________________________
std::unique_ptr<G4UIcommand> 
G4AnalysisMessengerHelper::CreateSetTitleCommand(G4UImessenger* messenger) const
{
  auto parId = new G4UIparameter("id", 'i', false);
  parId->SetGuidance(Update("OBJECT id"));
  parId->SetParameterRange("id>=0");

  auto parTitle = new G4UIparameter("title", 's', true);
  parTitle->SetGuidance(Update("OBJECT title"));
  parTitle->SetDefaultValue("none");

  std::unique_ptr<G4UIcommand> command(
    new G4UIcommand(Update("/analysis/HNTYPE_/setTitle"), messenger));
  command->SetGuidance(Update("Set title for the NDIM_D LOBJECT of given id"));
  command->SetParameter(parId);
  command->SetParameter(parTitle);
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}  


//_____________________________________________________________________________
std::unique_ptr<G4UIcommand> 
G4AnalysisMessengerHelper::CreateSetBinsCommand(const G4String& axis,
                                                G4UImessenger* messenger) const
{
  auto parId = new G4UIparameter("id", 'i', false);
  parId->SetGuidance(Update( "OBJECT id"));
  parId->SetParameterRange("id>=0");
  
  auto parNbins = new G4UIparameter("nbins", 'i', false);
  parNbins->SetGuidance("Number of bins");
  
  auto parValMin = new G4UIparameter("valMin", 'd', false);
  parValMin->SetGuidance("Minimum value, expressed in unit");
  
  auto parValMax = new G4UIparameter("valMax", 'd', false);
  parValMax->SetGuidance("Maximum value, expressed in unit");
  
  auto parValUnit = new G4UIparameter("valUnit", 's', true);
  parValUnit->SetGuidance("The unit applied to filled values and valMin, valMax");
  parValUnit->SetDefaultValue("none");

  auto parValFcn = new G4UIparameter("valFcn", 's', true);
  parValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnGuidance = "The function applied to filled values (log, log10, exp, none).\n";
  fcnGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnGuidance += "but none value should be used instead.";
  parValFcn->SetGuidance(fcnGuidance);
  parValFcn->SetDefaultValue("none");

  auto parValBinScheme = new G4UIparameter("valBinScheme", 's', true);
  parValBinScheme->SetParameterCandidates("linear log");
  G4String binSchemeGuidance = "The binning scheme (linear, log).\n";
  binSchemeGuidance 
    += "Note that the unit and fcn parameters cannot be omitted in this case,\n";
  binSchemeGuidance += "but none value should be used instead.";
  parValBinScheme->SetGuidance(binSchemeGuidance);
  parValBinScheme->SetDefaultValue("linear");
 
  auto commandName = Update("/analysis/HNTYPE_/setUAXIS", axis);
  std::unique_ptr<G4UIcommand> command(
    new G4UIcommand(Update("/analysis/HNTYPE_/setUAXIS", axis), messenger));
  command->SetGuidance(Update("Set parameters for the NDIM_D LOBJECT of given id:"));
  command->SetGuidance(
    Update("  nAXISbins; AXISvalMin; AXISvalMax; AXISunit; AXISfunction; AXISbinScheme", axis));
  command->SetParameter(parId);
  command->SetParameter(parNbins);
  command->SetParameter(parValMin);
  command->SetParameter(parValMax);
  command->SetParameter(parValUnit);
  command->SetParameter(parValFcn);
  command->SetParameter(parValBinScheme);
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}  

//_____________________________________________________________________________
 std::unique_ptr<G4UIcommand>  
 G4AnalysisMessengerHelper::CreateSetValuesCommand(const G4String& axis,
                                                   G4UImessenger* messenger) const
{
  auto parId = new G4UIparameter("id", 'i', false);
  parId->SetGuidance(Update("OBJECT id"));
  parId->SetParameterRange("id>=0");

  auto parValMin = new G4UIparameter("valMin", 'd', false);
  parValMin->SetGuidance(Update("Minimum AXIS-value expressed in unit", axis));
  
  auto parValMax = new G4UIparameter("valMax", 'd', false);
  parValMax->SetGuidance(Update("Maximum AXIS-value expressed in unit", axis));
  
  auto parValUnit = new G4UIparameter("valUnit", 's', true);
  parValUnit->SetGuidance("The unit applied to filled values and valMin, valMax");
  parValUnit->SetDefaultValue("none");

  auto parValFcn = new G4UIparameter("valFcn", 's', true);
  parValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnGuidance = "The function applied to filled values (log, log10, exp, none).\n";
  fcnGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnGuidance += "but none value should be used instead.";
  parValFcn->SetGuidance(fcnGuidance);
  parValFcn->SetDefaultValue("none");

  std::unique_ptr<G4UIcommand> command(
    new G4UIcommand(Update("/analysis/HNTYPE_/setUAXIS", axis), messenger));
  command->SetGuidance(Update("Set parameters for the NDIM_D LOBJECT of #id:"));
  command->SetGuidance(
    Update("  AXISvalMin; AXISvalMax; AXISunit; AXISfunction", axis));
  command->SetParameter(parId);
  command->SetParameter(parValMin);
  command->SetParameter(parValMax);
  command->SetParameter(parValUnit);
  command->SetParameter(parValFcn);
  command->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  return command;
}  

//_____________________________________________________________________________
std::unique_ptr<G4UIcommand> 
G4AnalysisMessengerHelper::CreateSetAxisCommand(const G4String& axis,
                                                G4UImessenger* messenger) const
{
  auto parId = new G4UIparameter("id", 'i', false);
  parId->SetGuidance(Update("OBJECT id"));
  parId->SetParameterRange("id>=0");

  auto parAxis = new G4UIparameter("axis", 's', true);
  parAxis->SetGuidance(Update("Histogram AXIS-axis title", axis));
  parAxis->SetDefaultValue("none");

  std::unique_ptr<G4UIcommand> command( 
    new G4UIcommand(Update("/analysis/HNTYPE_/setUAXISaxis", axis), messenger));
  command->SetGuidance(Update("Set AXIS-axis title for the NDIM_D LOBJECT of given id", axis));
  command->SetParameter(parId);
  command->SetParameter(parAxis);
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}

//_____________________________________________________________________________
void G4AnalysisMessengerHelper::GetBinData(BinData& data, 
                                           std::vector<G4String>& parameters, 
                                           G4int& counter) const
{
  data.fNbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
  data.fVmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
  data.fVmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
  data.fSunit = parameters[counter++];
  data.fSfcn = parameters[counter++];
  data.fSbinScheme = parameters[counter++];
}

//_____________________________________________________________________________
void G4AnalysisMessengerHelper::GetValueData(ValueData& data, 
                                             std::vector<G4String>& parameters, 
                                             G4int& counter) const
{
  data.fVmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
  data.fVmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
  data.fSunit = parameters[counter++];
  data.fSfcn = parameters[counter++];
}

//_____________________________________________________________________________
void G4AnalysisMessengerHelper::WarnAboutParameters(G4UIcommand* command, 
                                                    G4int nofParameters) const
{
  G4ExceptionDescription description;
  description 
    << "Got wrong number of \"" << command->GetCommandName() 
    << "\" parameters: " << nofParameters
    << " instead of " << command->GetParameterEntries() 
    << " expected" << G4endl;
  G4String methodName(Update("G4UHNTYPE_Messenger::SetNewValue"));  
  G4Exception(methodName,
              "Analysis_W013", JustWarning, description);
}

//_____________________________________________________________________________
void G4AnalysisMessengerHelper::WarnAboutSetCommands() const
{
  G4ExceptionDescription description;
  description 
    << "Command setX, setY, setZ must be called sucessively in this order. " << G4endl
    << "Command was ignored." << G4endl;
  G4String methodName(Update("G4UHNTYPE_Messenger::SetNewValue"));  
  G4Exception(methodName,
              "Analysis_W013", JustWarning, description);
}


