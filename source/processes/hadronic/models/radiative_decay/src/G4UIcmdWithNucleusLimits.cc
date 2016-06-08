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
#include "G4UIcmdWithNucleusLimits.hh"
#include "g4std/strstream"
////////////////////////////////////////////////////////////////////////////////
//
G4UIcmdWithNucleusLimits::G4UIcmdWithNucleusLimits
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * intParamAMin = new G4UIparameter('i');
  SetParameter(intParamAMin);
  G4UIparameter * intParamAMax = new G4UIparameter('i');
  SetParameter(intParamAMax);
  G4UIparameter * intParamZMin = new G4UIparameter('i');
  SetParameter(intParamZMin);
  G4UIparameter * intParamZMax = new G4UIparameter('i');
  SetParameter(intParamZMax);
}

////////////////////////////////////////////////////////////////////////////////
//
G4UIcmdWithNucleusLimits::~G4UIcmdWithNucleusLimits()
{
  ;
}
////////////////////////////////////////////////////////////////////////////////
//
G4NucleusLimits G4UIcmdWithNucleusLimits::
  GetNewNucleusLimitsValue(G4String paramString)
{
  G4int aMin;
  G4int aMax;
  G4int zMin;
  G4int zMax;
  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> aMin >> aMax >> zMin >> zMax;
  return G4NucleusLimits(aMin,aMax,zMin,zMax);
}
////////////////////////////////////////////////////////////////////////////////
//
G4String G4UIcmdWithNucleusLimits::ConvertToString
(G4NucleusLimits defLimits)
{
  char st[100];
  G4std::ostrstream os(st,100);
  os << defLimits.GetAMin() << " " << defLimits.GetAMax()
     << defLimits.GetZMin() << " " << defLimits.GetZMax()<< '\0';
  G4String vl = st;
  return vl;
}                         
////////////////////////////////////////////////////////////////////////////////
//
void G4UIcmdWithNucleusLimits::SetParameterName
(const char * theNameAMin,const char * theNameAMax,const char * theNameZMin,
const char * theNameZMax,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamAMin = GetParameter(0);
  theParamAMin->SetParameterName(theNameAMin);
  theParamAMin->SetOmittable(omittable);
  theParamAMin->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamAMax = GetParameter(1);
  theParamAMax->SetParameterName(theNameAMax);
  theParamAMax->SetOmittable(omittable);
  theParamAMax->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZMin = GetParameter(2);
  theParamZMin->SetParameterName(theNameZMin);
  theParamZMin->SetOmittable(omittable);
  theParamZMin->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZMax = GetParameter(3);
  theParamZMax->SetParameterName(theNameZMax);
  theParamZMax->SetOmittable(omittable);
  theParamZMax->SetCurrentAsDefault(currentAsDefault);
}
////////////////////////////////////////////////////////////////////////////////
//
void G4UIcmdWithNucleusLimits::SetDefaultValue(G4NucleusLimits defLimits)
{
  G4UIparameter * theParamAMin = GetParameter(0);
  theParamAMin->SetDefaultValue(defLimits.GetAMin());
  G4UIparameter * theParamAMax = GetParameter(1);
  theParamAMax->SetDefaultValue(defLimits.GetAMax());
  G4UIparameter * theParamZMin = GetParameter(2);
  theParamZMin->SetDefaultValue(defLimits.GetZMin());
  G4UIparameter * theParamZMax = GetParameter(3);
  theParamZMax->SetDefaultValue(defLimits.GetZMax());
}






