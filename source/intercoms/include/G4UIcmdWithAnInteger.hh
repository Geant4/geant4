// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithAnInteger.hh,v 1.1 1999-01-07 16:09:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWithAnInteger_H
#define G4UIcmdWithAnInteger_H 1

#include "G4UIcommand.hh"

class G4UIcmdWithAnInteger : public G4UIcommand
{
  public:
    G4UIcmdWithAnInteger
    (const char * theCommandPath,G4UImessenger * theMessenger);
    
    G4int GetNewIntValue(G4String paramString);
    G4String ConvertToString(G4int intValue);
    void SetParameterName(const char * theName,G4bool omittable,
                          G4bool currentAsDefault=false);
    void SetDefaultValue(G4int defVal);
};

#endif
