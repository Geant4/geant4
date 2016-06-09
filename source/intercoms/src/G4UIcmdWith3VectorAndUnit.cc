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
// $Id: G4UIcmdWith3VectorAndUnit.cc,v 1.7 2004/05/16 20:42:37 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include <strstream>

G4UIcmdWith3VectorAndUnit::G4UIcmdWith3VectorAndUnit
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * dblParamX = new G4UIparameter('d');
  SetParameter(dblParamX);
  G4UIparameter * dblParamY = new G4UIparameter('d');
  SetParameter(dblParamY);
  G4UIparameter * dblParamZ = new G4UIparameter('d');
  SetParameter(dblParamZ);
  G4UIparameter * untParam = new G4UIparameter('s');
  SetParameter(untParam);
  untParam->SetParameterName("Unit");
}

G4ThreeVector G4UIcmdWith3VectorAndUnit::GetNew3VectorValue(const char* paramString)
{
  return ConvertToDimensioned3Vector(paramString);
}

G4ThreeVector G4UIcmdWith3VectorAndUnit::GetNew3VectorRawValue(const char* paramString)
{
  G4double vx;
  G4double vy;
  G4double vz;
  char unts[30];
  std::istrstream is((char*)paramString);
  is >> vx >> vy >> vz >> unts;
  return G4ThreeVector(vx,vy,vz);
}

G4double G4UIcmdWith3VectorAndUnit::GetNewUnitValue(const char* paramString)
{
  G4double vx;
  G4double vy;
  G4double vz;
  char unts[30];
  std::istrstream is((char*)paramString);
  is >> vx >> vy >> vz >> unts;
  G4String unt = unts;
  return ValueOf(unt);
}

G4String G4UIcmdWith3VectorAndUnit::ConvertToStringWithBestUnit(G4ThreeVector vec)
{
  G4UIparameter* unitParam = GetParameter(3);
  G4String canList = unitParam->GetParameterCandidates();
  G4Tokenizer candidateTokenizer(canList);
  G4String aToken = candidateTokenizer();
  char strg[120];
  std::ostrstream os(strg,120);
  os << G4BestUnit(vec,CategoryOf(aToken)) << '\0';

  G4String st = strg;
  return st;
}

G4String G4UIcmdWith3VectorAndUnit::ConvertToStringWithDefaultUnit(G4ThreeVector vec)
{
  G4UIparameter* unitParam = GetParameter(3);
  G4String st;
  if(unitParam->IsOmittable())
  { st = ConvertToString(vec,unitParam->GetDefaultValue()); }
  else
  { st = ConvertToStringWithBestUnit(vec); }
  return st;
}

void G4UIcmdWith3VectorAndUnit::SetParameterName
(const char * theNameX,const char * theNameY,const char * theNameZ,
G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamX = GetParameter(0);
  theParamX->SetParameterName(theNameX);
  theParamX->SetOmittable(omittable);
  theParamX->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamY = GetParameter(1);
  theParamY->SetParameterName(theNameY);
  theParamY->SetOmittable(omittable);
  theParamY->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZ = GetParameter(2);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWith3VectorAndUnit::SetDefaultValue(G4ThreeVector vec)
{
  G4UIparameter * theParamX = GetParameter(0);
  theParamX->SetDefaultValue(vec.x());
  G4UIparameter * theParamY = GetParameter(1);
  theParamY->SetDefaultValue(vec.y());
  G4UIparameter * theParamZ = GetParameter(2);
  theParamZ->SetDefaultValue(vec.z());
}

void G4UIcmdWith3VectorAndUnit::SetUnitCategory(const char * unitCategory)
{
  SetUnitCandidates(UnitsList(unitCategory));
}

void G4UIcmdWith3VectorAndUnit::SetUnitCandidates(const char * candidateList)
{
  G4UIparameter * untParam = GetParameter(3);
  G4String canList = candidateList;
  untParam->SetParameterCandidates(canList);
}

void G4UIcmdWith3VectorAndUnit::SetDefaultUnit(const char * defUnit)
{
  G4UIparameter * untParam = GetParameter(3);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
  SetUnitCategory(CategoryOf(defUnit));
}

