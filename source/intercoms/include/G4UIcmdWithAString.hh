// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithAString.hh,v 1.1 1999-01-07 16:09:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWithAString_H
#define G4UIcmdWithAString_H 1

#include "G4UIcommand.hh"

class G4UIcmdWithAString : public G4UIcommand
{
  public:
    G4UIcmdWithAString
    (const char * theCommandPath,G4UImessenger * theMessenger);
    
    void SetParameterName(const char * theName,G4bool omittable,
                          G4bool currentAsDefault=false);
    void SetCandidates(const char * candidateList);
    void SetDefaultValue(const char * defVal);
};

#endif
