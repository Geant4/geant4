//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4UIcmdWith3VectorAndUnit.hh 67965 2013-03-13 09:35:29Z gcosmo $
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
    virtual G4int DoIt(G4String parameterList);
    static G4ThreeVector GetNew3VectorValue(const char* paramString);
    //  Convert string which represents three double values and a unit to
    // G4ThreeVector. Values are converted to the Geant4 internal unit.
    static G4ThreeVector GetNew3VectorRawValue(const char* paramString);
    //  Convert string which represents three double values and a unit to
    // G4ThreeVector. Values are NOT converted to the Geant4 internal unit
    // but just as the given string.
    static G4double GetNewUnitValue(const char* paramString);
    //  Convert the unit string to the value of the unit. "paramString"
    // must contain three double values AND a unit string.
    G4String ConvertToStringWithBestUnit(G4ThreeVector vec);
    //  Convert a 3 vector value to a string of digits and unit. Best unit is
    // chosen from the unit category of default unit (in case SetDefaultUnit()
    //  is defined) or category defined by SetUnitCategory().
    G4String ConvertToStringWithDefaultUnit(G4ThreeVector vec);
    //  Convert a 3 vector value to a string of digits and unit. Best unit is
    // chosen from the category defined by SetUnitCategory() in case default
    // unit is not defined.
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
