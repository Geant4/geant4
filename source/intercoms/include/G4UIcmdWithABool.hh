// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithABool.hh,v 2.1 1998/07/12 02:59:28 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
//

#ifndef G4UIcmdWithABool_H
#define G4UIcmdWithABool_H 1

#include "G4UIcommand.hh"

class G4UIcmdWithABool : public G4UIcommand
{
  public:
    G4UIcmdWithABool
    (const char * theCommandPath,G4UImessenger * theMessenger);
    
    G4bool GetNewBoolValue(G4String paramString);
    G4String ConvertToString(G4bool intValue);
    void SetParameterName(const char * theName,G4bool omittable,
                          G4bool currentAsDefault=false);
    void SetDefaultValue(G4bool defVal);
};

#endif
