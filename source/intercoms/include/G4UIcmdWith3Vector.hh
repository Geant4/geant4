// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWith3Vector.hh,v 2.1 1998/07/12 02:59:27 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
//

#ifndef G4UIcmdWith3Vector_H
#define G4UIcmdWith3Vector_H 1

#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"

class G4UIcmdWith3Vector : public G4UIcommand
{
  public:
    G4UIcmdWith3Vector
    (const char * theCommandPath,G4UImessenger * theMessenger);
    
    G4ThreeVector GetNew3VectorValue(G4String paramString);
    G4String ConvertToString(G4ThreeVector vec);
    void SetParameterName(const char * theNameX,const char * theNameY,
      const char * theNameZ,G4bool omittable,G4bool currentAsDefault=false);
    void SetDefaultValue(G4ThreeVector defVal);
};

#endif
