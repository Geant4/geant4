// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithAnInteger.cc,v 1.1 1999-01-07 16:09:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4UIcmdWithAnInteger.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4UIcmdWithAnInteger::G4UIcmdWithAnInteger
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * intParam = new G4UIparameter('i');
  SetParameter(intParam);
}

G4int G4UIcmdWithAnInteger::GetNewIntValue(G4String paramString)
{
  G4int vl;
  const char* t = paramString;
  istrstream is((char*)t);
  is >> vl;
  return vl;
}

G4String G4UIcmdWithAnInteger::ConvertToString(G4int intValue)
{
  char st[20];
  ostrstream os(st,20);
  os << intValue << '\0';
  G4String vl = st;
  return vl;
}

void G4UIcmdWithAnInteger::SetParameterName
(const char * theName,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWithAnInteger::SetDefaultValue(G4int defVal)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}

