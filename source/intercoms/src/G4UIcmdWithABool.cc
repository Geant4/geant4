// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithABool.cc,v 1.1 1999-01-07 16:09:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4UIcmdWithABool.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4UIcmdWithABool::G4UIcmdWithABool
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * blParam = new G4UIparameter('b');
  SetParameter(blParam);
}

G4bool G4UIcmdWithABool::GetNewBoolValue(G4String paramString)
{
  G4String v = paramString;
  v.toUpper();
  G4bool vl = false;
  if( v=="Y" || v=="YES" || v=="1" || v=="T" || v=="TRUE" )
  { vl = true; }
  return vl;
}

G4String G4UIcmdWithABool::ConvertToString(G4bool blValue)
{
  G4String vl = "0";
  if(blValue) vl = "1";
  return vl;
}

void G4UIcmdWithABool::SetParameterName
(const char * theName,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWithABool::SetDefaultValue(G4bool defVal)
{
  G4UIparameter * theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}

