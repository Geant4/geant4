//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4UIcmdWithADoubleAndUnit.hh,v 1.5 2002-04-26 22:03:35 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef G4UIcmdWithADoubleAndUnit_H
#define G4UIcmdWithADoubleAndUnit_H 1

#include "G4UIcommand.hh"

// class description:
//  A concrete class of G4UIcommand. The command defined by this class
// takes a double value and a unit string.
//  General information of G4UIcommand is given in G4UIcommand.hh.

class G4UIcmdWithADoubleAndUnit : public G4UIcommand
{
  public: // with description
    G4UIcmdWithADoubleAndUnit
    (const char * theCommandPath,G4UImessenger * theMessenger);
    //  Constructor. The command string with full path directory
    // and the pointer to the messenger must be given.    
    G4double GetNewDoubleValue(const char* paramString);
    //  Convert string which represents a double value and a unit to
    // double. Value is converted to the Geant4 internal unit.
    G4double GetNewDoubleRawValue(const char* paramString);
    //  Convert string which represents a double value and a unit to
    // double. Value is NOT converted to the Geant4 internal unit
    // but just as the given string.
    G4double GetNewUnitValue(const char* paramString);
    //  Convert the unit string to the value of the unit. "paramString"
    // must contain a double value AND a unit string.
    G4String ConvertToString(G4double dblValue,const char * unitName);
    //  Convert a double value and the unit to a string which represents a
    // double value and the unit. This method must be used by the messenger
    // for its GetCurrentValues() method. Value of the double will be
    // devided by the value of the unit.
    void SetParameterName(const char * theName,G4bool omittable,
                          G4bool currentAsDefault=false);
    //  Set the parameter name for double parameterxs. Name is used by
    // the range checking routine.
    //  If "omittable" is set as true, the user of this command can ommit
    // the value when he/she applies the command. If "omittable" is false,
    // the user must supply a value.
    //  "currentAsDefault" flag is valid only if "omittable" is true. If this
    // flag is true, the current value is used as the default value when the 
    // user ommit the double parameter. If this flag is false, the value
    // given by the next SetDefaultValue() method is used.
    void SetDefaultValue(G4double defVal);
    //  Set the default value of the parameter. This default value is used
    // when the user of this command ommits the parameter value, and
    // "ommitable" is true and "curreutAsDefault" is false.
    void SetUnitCategory(const char * unitCategory);
    void SetUnitCandidates(const char * candidateList);
    void SetDefaultUnit(const char * defUnit);
    //  These three methods must be used alternatively.
    //  The user cannot ommit the unit as the second parameter of the command if
    // SetUnitCategory() or SetUnitCandidates() is used, while the unit defined
    // by SetDefaultUnit() method is used as the default unit so that the user can
    // ommits the second parameter.
    //  SetUnitCategory() defines the category of the units which will be accepted.
    // The available categories can be found in G4SystemOfUnits.hh in global category.
    // Only the units categorized in the given category are accepted as the second
    // parameter of the command.
    //  SetUnitCandidates() defines the candidates of units. Units listed in the
    // argument of this method must be separated by space(s). Only the units listed
    // in the candidate list are accepted as the second parameter of the command.
    //  SetDefaultUnit() defines the default unit and also it defines the category
    // of the allowed units. Thus only the units categorized as the given default
    // unit will be accepted.
};

#endif
