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
// G4UIcmdWithADoubleAndUnit
//
// Author: M.Asai, 1998
// --------------------------------------------------------------------

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommandStatus.hh"
#include <sstream>
#include <vector>

// --------------------------------------------------------------------
G4UIcmdWithADoubleAndUnit::G4UIcmdWithADoubleAndUnit(
  const char* theCommandPath, G4UImessenger* theMessenger)
  : G4UIcommand(theCommandPath, theMessenger)
{
  auto* dblParam = new G4UIparameter('d');
  SetParameter(dblParam);
  auto* untParam = new G4UIparameter('s');
  untParam->SetParameterName("Unit");
  SetParameter(untParam);
  SetCommandType(WithADoubleAndUnitCmd);
}

// --------------------------------------------------------------------
G4int G4UIcmdWithADoubleAndUnit::DoIt(G4String parameterList)
{
  std::vector<G4String> token_vector;
  G4Tokenizer tkn(parameterList);
  G4String str;
  while(!(str = tkn()).empty())
  {
    token_vector.push_back(str);
  }

  // convert a value in default unit
  G4String converted_parameter;
  G4String default_unit = GetParameter(1)->GetDefaultValue();
  if(!default_unit.empty() && token_vector.size() >= 2)
  {
    if(CategoryOf(token_vector[1]) != CategoryOf(default_unit))
    {
      return fParameterOutOfCandidates + 1;
    }
    G4double value_given   = ValueOf(token_vector[1]);
    G4double value_default = ValueOf(default_unit);
    G4double value =
      ConvertToDouble(token_vector[0]) * value_given / value_default;
    // reconstruct parameter list
    converted_parameter += ConvertToString(value);
    converted_parameter += " ";
    converted_parameter += default_unit;
    for(std::size_t i = 2; i < token_vector.size(); ++i)
    {
      converted_parameter += " ";
      converted_parameter += token_vector[i];
    }
  }
  else
  {
    converted_parameter = parameterList;
  }

  return G4UIcommand::DoIt(converted_parameter);
}

// --------------------------------------------------------------------
G4double G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(const char* paramString)
{
  return ConvertToDimensionedDouble(paramString);
}

// --------------------------------------------------------------------
G4double G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(
  const char* paramString)
{
  G4double vl;
  char unts[30];

  std::istringstream is(paramString);
  is >> vl >> unts;

  return vl;
}

// --------------------------------------------------------------------
G4double G4UIcmdWithADoubleAndUnit::GetNewUnitValue(const char* paramString)
{
  G4double vl;
  char unts[30];

  std::istringstream is(paramString);
  is >> vl >> unts;
  G4String unt = unts;

  return ValueOf(unt);
}

// --------------------------------------------------------------------
G4String G4UIcmdWithADoubleAndUnit::ConvertToStringWithBestUnit(G4double val)
{
  G4UIparameter* unitParam = GetParameter(1);
  G4String canList         = unitParam->GetParameterCandidates();
  G4Tokenizer candidateTokenizer(canList);
  G4String aToken = candidateTokenizer();
  std::ostringstream os;
  os << G4BestUnit(val, CategoryOf(aToken));

  G4String st = os.str();
  return st;
}

// --------------------------------------------------------------------
G4String G4UIcmdWithADoubleAndUnit::ConvertToStringWithDefaultUnit(G4double val)
{
  G4UIparameter* unitParam = GetParameter(1);
  G4String st;
  if(unitParam->IsOmittable())
  {
    st = ConvertToString(val, unitParam->GetDefaultValue());
  }
  else
  {
    st = ConvertToStringWithBestUnit(val);
  }
  return st;
}

// --------------------------------------------------------------------
void G4UIcmdWithADoubleAndUnit::SetParameterName(const char* theName,
                                                 G4bool omittable,
                                                 G4bool currentAsDefault)
{
  G4UIparameter* theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

// --------------------------------------------------------------------
void G4UIcmdWithADoubleAndUnit::SetDefaultValue(G4double defVal)
{
  G4UIparameter* theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}

// --------------------------------------------------------------------
void G4UIcmdWithADoubleAndUnit::SetUnitCategory(const char* unitCategory)
{
  SetUnitCandidates(UnitsList(unitCategory));
}

// --------------------------------------------------------------------
void G4UIcmdWithADoubleAndUnit::SetUnitCandidates(const char* candidateList)
{
  G4UIparameter* untParam = GetParameter(1);
  G4String canList        = candidateList;
  untParam->SetParameterCandidates(canList);
}

// --------------------------------------------------------------------
void G4UIcmdWithADoubleAndUnit::SetDefaultUnit(const char* defUnit)
{
  G4UIparameter* untParam = GetParameter(1);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
  SetUnitCategory(CategoryOf(defUnit));
}
