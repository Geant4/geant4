// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWithABool.hh,v 1.3 1999-12-15 14:50:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWithABool_H
#define G4UIcmdWithABool_H 1

#include "G4UIcommand.hh"

// class description:
//  A concrete class of G4UIcommand. The command defined by this class
// takes a boolean value. Boolean value can be the following notations.
//    TRUE : 
//       1 t T true TRUE
//    FALSE : 
//       0 f F false FALSE
//  General information of G4UIcommand is given in G4UIcommand.hh.

class G4UIcmdWithABool : public G4UIcommand
{
  public: // with description
    G4UIcmdWithABool
    (const char * theCommandPath,G4UImessenger * theMessenger);
    //  Constructor. The command string with full path directory
    // and the pointer to the messenger must be given.
    G4bool GetNewBoolValue(G4String paramString);
    //  Convert string which represents a boolean value to G4bool.
    G4String ConvertToString(G4bool intValue);
    //  Convert a boolean value to a string. This method must be used by 
    // the messenger for its GetCurrentValues() method.
    void SetParameterName(const char * theName,G4bool omittable,
                          G4bool currentAsDefault=false);
    //  Set the parameter name for a boolean parameter.
    //  If "omittable" is set as true, the user of this command can ommit
    // the value when he/she applies the command. If "omittable" is false,
    // the user must supply a boolean value.
    //  "currentAsDefault" flag is valid only if "omittable" is true. If this
    // flag is true, the current value is used as the default value when the 
    // user ommit the parameter. If this flag is false, the value given by the 
    // next SetDefaultValue() method is used.
    void SetDefaultValue(G4bool defVal);
    //  Set the default value of the parameter. This default value is used
    // when the user of this command ommits the parameter value, and
    // "ommitable" is true and "currentAsDefault" is false.
};

#endif
