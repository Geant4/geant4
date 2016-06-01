// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithADoubleAndUnit.hh,v 2.1 1998/07/12 02:59:30 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
//

#ifndef G4UIcmdWithADoubleAndUnit_H
#define G4UIcmdWithADoubleAndUnit_H 1

#include "G4UIcommand.hh"

class G4UIcmdWithADoubleAndUnit : public G4UIcommand
{
  public:
    G4UIcmdWithADoubleAndUnit
    (const char * theCommandPath,G4UImessenger * theMessenger);
    
    G4double GetNewDoubleValue(G4String paramString);
    G4double GetNewDoubleRawValue(G4String paramString);
    G4double GetNewUnitValue(G4String paramString);
    G4String ConvertToString(G4double dblValue,const char * unitName);
    void SetParameterName(const char * theName,G4bool omittable,
                          G4bool currentAsDefault=false);
    void SetDefaultValue(G4double defVal);
    void SetUnitCategory(const char * unitCategory);
    void SetUnitCandidates(const char * candidateList);
    void SetDefaultUnit(const char * defUnit);
};

#endif
