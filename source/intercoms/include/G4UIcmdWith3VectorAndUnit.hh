// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcmdWith3VectorAndUnit.hh,v 1.2 1999-10-29 06:06:43 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWith3VectorAndUnit_H
#define G4UIcmdWith3VectorAndUnit_H 1

#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"

// class description:
//  A concrete class of G4UIcommand. The command defined by this class
// takes three double values and a unit.
//  General information of G4UIcommand is given in G4UIcommand.hh.

class G4UIcmdWith3VectorAndUnit : public G4UIcommand
{
  public: // with description
    G4UIcmdWith3VectorAndUnit
    (const char * theCommandPath,G4UImessenger * theMessenger);
    //  Constructor. The command string with full path directory
    // and the pointer to the messenger must be given.
    G4ThreeVector GetNew3VectorValue(G4String paramString);
    //  Convert string which represents three double values and a unit to
    // G4ThreeVector. Values are converted to the Geant4 internal unit.
    G4ThreeVector GetNew3VectorRawValue(G4String paramString);
    //  Convert string which represents three double values and a unit to
    // G4ThreeVector. Values are NOT converted to the Geant4 internal unit
    // but just as the given string.
    G4double GetNewUnitValue(G4String paramString);
    //  Convert the unit string to the value of the unit. "paramString"
    // must contain three double values AND a unit string.
    G4String ConvertToString(G4ThreeVector vec,const char * unitName);
    //  Convert G4ThreeVector and the unit to a string which represents three
    // double values and the unit. This method must be used by the messenger
    // for its GetCurrentValues() method. Values of the three vector will be
    // devided by the value of the unit.
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
    void SetUnitCategory(const char * unitCategory);
    void SetUnitCandidates(const char * candidateList);
    void SetDefaultUnit(const char * defUnit);
    //  These three methods must be used alternatively.
    //  The user cannot ommit the unit as the fourth parameter of the command if
    // SetUnitCategory() or SetUnitCandidates() is used, while the unit defined
    // by SetDefaultUnit() method is used as the default unit so that the user can
    // ommits the fourth parameter.
    //  SetUnitCategory() defines the category of the units which will be accepted.
    // The available categories can be found in G4SystemOfUnits.hh in global category.
    // Only the units categorized in the given category are accepted as the fourth
    // parameter of the command.
    //  SetUnitCandidates() defines the candidates of units. Units listed in the
    // argument of this method must be separated by space(s). Only the units listed
    // in the candidate list are accepted as the fourth parameter of the command.
    //  SetDefaultUnit() defines the default unit and also it defines the category
    // of the allowed units. Thus only the units categorized as the given default
    // unit will be accepted.
};

#endif
