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
// G4UIcmdWithALongInt
//
// Author: M.Asai, 2020
// --------------------------------------------------------------------

#include "G4UIcmdWithALongInt.hh"

// --------------------------------------------------------------------
G4UIcmdWithALongInt::G4UIcmdWithALongInt(const char* theCommandPath,
                                           G4UImessenger* theMessenger)
  : G4UIcommand(theCommandPath, theMessenger)
{
  auto* longParam = new G4UIparameter('l');
  SetParameter(longParam);
  SetCommandType(WithALongIntCmd);
}

// --------------------------------------------------------------------
G4long G4UIcmdWithALongInt::GetNewLongIntValue(const char* paramString)
{
  return ConvertToLongInt(paramString);
}

// --------------------------------------------------------------------
void G4UIcmdWithALongInt::SetParameterName(const char* theName,
                                            G4bool omittable,
                                            G4bool currentAsDefault)
{
  G4UIparameter* theParam = GetParameter(0);
  theParam->SetParameterName(theName);
  theParam->SetOmittable(omittable);
  theParam->SetCurrentAsDefault(currentAsDefault);
}

// --------------------------------------------------------------------
void G4UIcmdWithALongInt::SetDefaultValue(G4long defVal)
{
  G4UIparameter* theParam = GetParameter(0);
  theParam->SetDefaultValue(defVal);
}
