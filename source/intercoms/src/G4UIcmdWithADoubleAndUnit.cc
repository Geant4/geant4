// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithADoubleAndUnit.cc,v 1.1 1999-01-07 16:09:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4UIcmdWithADoubleAndUnit.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

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

G4double G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(G4String paramString)
{
  G4double vl;
  char unts[30];
  
  const char* t = paramString;
  istrstream is((char*)t);
  is >> vl >> unts;
  G4String unt = unts;
  
  return (vl*ValueOf(unt));
}

G4double G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(G4String paramString)
{
  G4double vl;
  char unts[30];
  
  const char* t = paramString;
  istrstream is((char*)t);
  is >> vl >> unts;
  
  return vl;
}

G4double G4UIcmdWithADoubleAndUnit::GetNewUnitValue(G4String paramString)
{
  G4double vl;
  char unts[30];
  
  const char* t = paramString;
  istrstream is((char*)t);
  is >> vl >> unts;
  G4String unt = unts;
  
  return ValueOf(unt);
}

G4String G4UIcmdWithADoubleAndUnit::ConvertToString
(G4double dblValue,const char * unitName)
{
  G4String unt = unitName;
  G4double uv = ValueOf(unitName);
  
  char st[50];
  ostrstream os(st,50);
  os << dblValue/uv << " " << unitName << '\0';
  G4String vl = st;
  return vl;
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

