// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithADouble.cc,v 1.1 1999-01-07 16:09:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4UIcmdWithADouble.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4UIcmdWithADouble::G4UIcmdWithADouble
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * dblParam = new G4UIparameter('d');
  SetParameter(dblParam);
}

G4double G4UIcmdWithADouble::GetNewDoubleValue(G4String paramString)
{
  G4double vl;
  const char* t = paramString;
  istrstream is((char*)t);
  is >> vl;
  return vl;
}

G4String G4UIcmdWithADouble::ConvertToString(G4double dblValue)
{
  char st[20];
  ostrstream os(st,20);
  os << dblValue << '\0';
  G4String vl = st;
  return vl;
}

void G4UIcmdWithADouble::SetParameterName
(const char * theName,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWithADouble::SetDefaultValue(G4double defVal)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}

