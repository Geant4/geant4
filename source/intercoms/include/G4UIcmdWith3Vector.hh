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
// G4UIcmdWith3Vector
//
// Class description:
//
// A concrete class of G4UIcommand. The command defined by this class
// takes three double values.
// General information of G4UIcommand is given in G4UIcommand.hh

// Author: M.Asai, 1998
// --------------------------------------------------------------------
#ifndef G4UIcmdWith3Vector_hh
#define G4UIcmdWith3Vector_hh 1

#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"

class G4UIcmdWith3Vector : public G4UIcommand
{
  public:

    G4UIcmdWith3Vector(const char* theCommandPath, G4UImessenger* theMessenger);
      // Constructor. The command string with full path directory
      // and the pointer to the messenger must be given

    static G4ThreeVector GetNew3VectorValue(const char* paramString);
      // Convert string which represents three double values to
      // G4ThreeVector

    void SetParameterName(const char* theNameX, const char* theNameY,
                          const char* theNameZ, G4bool omittable,
                          G4bool currentAsDefault = false);
      // Set the parameter names for three parameters. Names are used by
      // the range checking function.
      // If "omittable" is set as true, the user of this command can omit
      // the value(s) when the command is applied. If "omittable" is false,
      // the user must supply all three values.
      // "currentAsDefault" flag is valid only if "omittable" is true. If this
      // flag is true, the current values are used as the default values when
      // the user omits some of the parameters. If this flag is false, the
      // values given by the next SetDefaultValue() method are used

    void SetDefaultValue(const G4ThreeVector& defVal);
    // Set the default values of the parameters. These default values are
    // used when the user of this command omits some of the parameter values,
    // and "omittable" is true and "currentAsDefault" is false
};

#endif
