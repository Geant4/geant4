// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithAString.cc,v 1.2 1999-12-15 14:50:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4UIcmdWithAString.hh"
#include "g4std/strstream"

G4UIcmdWithAString::G4UIcmdWithAString
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * strParam = new G4UIparameter('s');
  SetParameter(strParam);
}

void G4UIcmdWithAString::SetParameterName
(const char * theName,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWithAString::SetCandidates(const char * candidateList)
{
  G4UIparameter * theParam = GetParameter(0);
  G4String canList = candidateList;
  theParam->SetParameterCandidates(canList);
}

void G4UIcmdWithAString::SetDefaultValue(const char * defVal)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}

