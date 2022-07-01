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
// G4UIcmdWith3VectorAndUnit
//
// Author: M.Asai, 1998
// --------------------------------------------------------------------

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommandStatus.hh"
#include <sstream>

// --------------------------------------------------------------------
G4UIcmdWith3VectorAndUnit::G4UIcmdWith3VectorAndUnit(
  const char* theCommandPath, G4UImessenger* theMessenger)
  : G4UIcommand(theCommandPath, theMessenger)
{
  auto* dblParamX = new G4UIparameter('d');
  SetParameter(dblParamX);
  auto* dblParamY = new G4UIparameter('d');
  SetParameter(dblParamY);
  auto* dblParamZ = new G4UIparameter('d');
  SetParameter(dblParamZ);
  auto* untParam = new G4UIparameter('s');
  untParam->SetParameterName("Unit");
  SetParameter(untParam);
  SetCommandType(With3VectorAndUnitCmd);
}

// --------------------------------------------------------------------
G4int G4UIcmdWith3VectorAndUnit::DoIt(G4String parameterList)
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
  G4String default_unit = GetParameter(3)->GetDefaultValue();
  if(!default_unit.empty() && token_vector.size() >= 4)
  {
    if(CategoryOf(token_vector[3]) != CategoryOf(default_unit))
    {
      return fParameterOutOfCandidates + 3;
    }
    G4double value_given   = ValueOf(token_vector[3]);
    G4double value_default = ValueOf(default_unit);
    G4double x = ConvertToDouble(token_vector[0]) * value_given / value_default;
    G4double y = ConvertToDouble(token_vector[1]) * value_given / value_default;
    G4double z = ConvertToDouble(token_vector[2]) * value_given / value_default;

    // reconstruct parameter list
    converted_parameter += ConvertToString(x);
    converted_parameter += " ";
    converted_parameter += ConvertToString(y);
    converted_parameter += " ";
    converted_parameter += ConvertToString(z);
    converted_parameter += " ";
    converted_parameter += default_unit;
    for(std::size_t i = 4; i < token_vector.size(); ++i)
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
G4ThreeVector G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(
  const char* paramString)
{
  return ConvertToDimensioned3Vector(paramString);
}

// --------------------------------------------------------------------
G4ThreeVector G4UIcmdWith3VectorAndUnit::GetNew3VectorRawValue(
  const char* paramString)
{
  G4double vx;
  G4double vy;
  G4double vz;
  char unts[30];
  std::istringstream is(paramString);
  is >> vx >> vy >> vz >> unts;
  return G4ThreeVector(vx, vy, vz);
}

// --------------------------------------------------------------------
G4double G4UIcmdWith3VectorAndUnit::GetNewUnitValue(const char* paramString)
{
  G4double vx;
  G4double vy;
  G4double vz;
  char unts[30];
  std::istringstream is(paramString);
  is >> vx >> vy >> vz >> unts;
  G4String unt = unts;
  return ValueOf(unt);
}

// --------------------------------------------------------------------
G4String G4UIcmdWith3VectorAndUnit::ConvertToStringWithBestUnit(
  const G4ThreeVector& vec)
{
  G4UIparameter* unitParam = GetParameter(3);
  G4String canList         = unitParam->GetParameterCandidates();
  G4Tokenizer candidateTokenizer(canList);
  G4String aToken = candidateTokenizer();

  std::ostringstream os;
  os << G4BestUnit(vec, CategoryOf(aToken));
  G4String st = os.str();

  return st;
}

// --------------------------------------------------------------------
G4String G4UIcmdWith3VectorAndUnit::ConvertToStringWithDefaultUnit(
  const G4ThreeVector& vec)
{
  G4UIparameter* unitParam = GetParameter(3);
  G4String st;
  if(unitParam->IsOmittable())
  {
    st = ConvertToString(vec, unitParam->GetDefaultValue());
  }
  else
  {
    st = ConvertToStringWithBestUnit(vec);
  }
  return st;
}

// --------------------------------------------------------------------
void G4UIcmdWith3VectorAndUnit::SetParameterName(const char* theNameX,
                                                 const char* theNameY,
                                                 const char* theNameZ,
                                                 G4bool omittable,
                                                 G4bool currentAsDefault)
{
  G4UIparameter* theParamX = GetParameter(0);
  theParamX->SetParameterName(theNameX);
  theParamX->SetOmittable(omittable);
  theParamX->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter* theParamY = GetParameter(1);
  theParamY->SetParameterName(theNameY);
  theParamY->SetOmittable(omittable);
  theParamY->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter* theParamZ = GetParameter(2);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
}

// --------------------------------------------------------------------
void G4UIcmdWith3VectorAndUnit::SetDefaultValue(const G4ThreeVector& vec)
{
  G4UIparameter* theParamX = GetParameter(0);
  theParamX->SetDefaultValue(vec.x());
  G4UIparameter* theParamY = GetParameter(1);
  theParamY->SetDefaultValue(vec.y());
  G4UIparameter* theParamZ = GetParameter(2);
  theParamZ->SetDefaultValue(vec.z());
}

// --------------------------------------------------------------------
void G4UIcmdWith3VectorAndUnit::SetUnitCategory(const char* unitCategory)
{
  SetUnitCandidates(UnitsList(unitCategory));
}

// --------------------------------------------------------------------
void G4UIcmdWith3VectorAndUnit::SetUnitCandidates(const char* candidateList)
{
  G4UIparameter* untParam = GetParameter(3);
  G4String canList        = candidateList;
  untParam->SetParameterCandidates(canList);
}

// --------------------------------------------------------------------
void G4UIcmdWith3VectorAndUnit::SetDefaultUnit(const char* defUnit)
{
  G4UIparameter* untParam = GetParameter(3);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
  SetUnitCategory(CategoryOf(defUnit));
}
