// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWith3Vector.hh,v 1.2 1999-10-29 06:06:43 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWith3Vector_H
#define G4UIcmdWith3Vector_H 1

#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"

// class description:
//  A concrete class of G4UIcommand. The command defined by this class
// takes three double values.
//  General information of G4UIcommand is given in G4UIcommand.hh.

class G4UIcmdWith3Vector : public G4UIcommand
{
  public: // with description
    G4UIcmdWith3Vector
    (const char * theCommandPath,G4UImessenger * theMessenger);
    //  Constructor. The command string with full path directory
    // and the pointer to the messenger must be given.
    G4ThreeVector GetNew3VectorValue(G4String paramString);
    //  Convert string which represents three double values to
    // G4ThreeVector.
    G4String ConvertToString(G4ThreeVector vec);
    //  Convert G4ThreeVector to a string which represents three
    // double values. This method must be used by the messenger
    // for its GetCurrentValues() method.
    void SetParameterName(const char * theNameX,const char * theNameY,
      const char * theNameZ,G4bool omittable,G4bool currentAsDefault=false);
    //  Set the parameter names for three parameters. Names are used by
    // the range checking routine.
    //  If "omittable" is set as true, the user of this command can ommit
    // the value(s) when he/she applies the command. If "omittable" is false,
    // the user must supply all three values.
    //  "currentAsDefault" flag is valid only if "omittable" is true. If this
    // flag is true, the current values are used as the default values when the 
    // user ommit some of the parameters. If this flag is false, the values
    // given by the next SetDefaultValue() method are used. 
    void SetDefaultValue(G4ThreeVector defVal);
    //  Set the default values of the parameters. These default values are used
    // when the user of this command ommits some of the parameter values, and
    // "ommitable" is true and "currentAsDefault" is false.
};

#endif
