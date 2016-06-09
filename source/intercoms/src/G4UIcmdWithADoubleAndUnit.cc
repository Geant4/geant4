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
// $Id: G4UIcmdWithADoubleAndUnit.cc,v 1.7 2004/05/16 20:42:37 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include <strstream>

G4UIcmdWithADoubleAndUnit::G4UIcmdWithADoubleAndUnit
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * dblParam = new G4UIparameter('d');
  SetParameter(dblParam);
  G4UIparameter * untParam = new G4UIparameter('s');
  SetParameter(untParam);
  untParam->SetParameterName("Unit");
}

G4double G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(const char* paramString)
{
  return ConvertToDimensionedDouble(paramString);
}

G4double G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(const char* paramString)
{
  G4double vl;
  char unts[30];
  
  std::istrstream is((char*)paramString);
  is >> vl >> unts;
  
  return vl;
}

G4double G4UIcmdWithADoubleAndUnit::GetNewUnitValue(const char* paramString)
{
  G4double vl;
  char unts[30];
  
  std::istrstream is((char*)paramString);
  is >> vl >> unts;
  G4String unt = unts;
  
  return ValueOf(unt);
}

G4String G4UIcmdWithADoubleAndUnit::ConvertToStringWithBestUnit(G4double val)
{
  G4UIparameter* unitParam = GetParameter(1);
  G4String canList = unitParam->GetParameterCandidates();
  G4Tokenizer candidateTokenizer(canList);
  G4String aToken = candidateTokenizer();
  char strg[60];
  std::ostrstream os(strg,60);
  os << G4BestUnit(val,CategoryOf(aToken)) << '\0';

  G4String st = strg;
  return st;
}

G4String G4UIcmdWithADoubleAndUnit::ConvertToStringWithDefaultUnit(G4double val)
{
  G4UIparameter* unitParam = GetParameter(1);
  G4String st;
  if(unitParam->IsOmittable())
  { st = ConvertToString(val,unitParam->GetDefaultValue()); }
  else
  { st = ConvertToStringWithBestUnit(val); }
  return st;
}

void G4UIcmdWithADoubleAndUnit::SetParameterName
(const char * theName,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWithADoubleAndUnit::SetDefaultValue(G4double defVal)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}

void G4UIcmdWithADoubleAndUnit::SetUnitCategory(const char * unitCategory)
{
  SetUnitCandidates(UnitsList(unitCategory));
}

void G4UIcmdWithADoubleAndUnit::SetUnitCandidates(const char * candidateList)
{
  G4UIparameter * untParam = GetParameter(1);
  G4String canList = candidateList;
  untParam->SetParameterCandidates(canList);
}

void G4UIcmdWithADoubleAndUnit::SetDefaultUnit(const char * defUnit)
{
  G4UIparameter * untParam = GetParameter(1);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
  SetUnitCategory(CategoryOf(defUnit));
}

