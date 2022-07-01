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
// Author: M.Asai, 1998
// --------------------------------------------------------------------

#include "G4UIcmdWith3Vector.hh"

// --------------------------------------------------------------------
G4UIcmdWith3Vector::G4UIcmdWith3Vector(const char* theCommandPath,
                                       G4UImessenger* theMessenger)
  : G4UIcommand(theCommandPath, theMessenger)
{
  auto* dblParamX = new G4UIparameter('d');
  SetParameter(dblParamX);
  auto* dblParamY = new G4UIparameter('d');
  SetParameter(dblParamY);
  auto* dblParamZ = new G4UIparameter('d');
  SetParameter(dblParamZ);
  SetCommandType(With3VectorCmd);
}

// --------------------------------------------------------------------
G4ThreeVector G4UIcmdWith3Vector::GetNew3VectorValue(const char* paramString)
{
  return ConvertTo3Vector(paramString);
}

// --------------------------------------------------------------------
void G4UIcmdWith3Vector::SetParameterName(const char* theNameX,
                                          const char* theNameY,
                                          const char* theNameZ,
                                          G4bool omittable,
                                          G4bool currentAsDefault)
{
  G4UIparameter* theParamX = GetParameter(0);
  theParamX->SetParameterName(theNameX);
  theParamX->SetOmittable(omittable);
  theParamX->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter* theParamY = GetParameter(1);
  theParamY->SetParameterName(theNameY);
  theParamY->SetOmittable(omittable);
  theParamY->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter* theParamZ = GetParameter(2);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
}

// --------------------------------------------------------------------
void G4UIcmdWith3Vector::SetDefaultValue(const G4ThreeVector& vec)
{
  G4UIparameter* theParamX = GetParameter(0);
  theParamX->SetDefaultValue(vec.x());
  G4UIparameter* theParamY = GetParameter(1);
  theParamY->SetDefaultValue(vec.y());
  G4UIparameter* theParamZ = GetParameter(2);
  theParamZ->SetDefaultValue(vec.z());
}
