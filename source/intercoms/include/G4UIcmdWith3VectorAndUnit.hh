// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWith3VectorAndUnit.hh,v 1.1 1999-01-07 16:09:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWith3VectorAndUnit_H
#define G4UIcmdWith3VectorAndUnit_H 1

#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"

class G4UIcmdWith3VectorAndUnit : public G4UIcommand
{
  public:
    G4UIcmdWith3VectorAndUnit
    (const char * theCommandPath,G4UImessenger * theMessenger);
    
    G4ThreeVector GetNew3VectorValue(G4String paramString);
    G4ThreeVector GetNew3VectorRawValue(G4String paramString);
    G4double GetNewUnitValue(G4String paramString);
    G4String ConvertToString(G4ThreeVector vec,const char * unitName);
    void SetParameterName(const char * theNameX,const char * theNameY,
      const char * theNameZ,G4bool omittable,G4bool currentAsDefault=false);
    void SetDefaultValue(G4ThreeVector defVal);
    void SetUnitCategory(const char * unitCategory);
    void SetUnitCandidates(const char * candidateList);
    void SetDefaultUnit(const char * defUnit);
};

#endif
