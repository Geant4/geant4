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
// $Id: G4UIcmdWith3Vector.cc,v 1.7 2005/10/26 06:10:22 kmura Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//

#include "G4UIcmdWith3Vector.hh"

G4UIcmdWith3Vector::G4UIcmdWith3Vector
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * dblParamX = new G4UIparameter('d');
  SetParameter(dblParamX);
  G4UIparameter * dblParamY = new G4UIparameter('d');
  SetParameter(dblParamY);
  G4UIparameter * dblParamZ = new G4UIparameter('d');
  SetParameter(dblParamZ);
}

G4ThreeVector G4UIcmdWith3Vector::GetNew3VectorValue(const char* paramString)
{
  return ConvertTo3Vector(paramString);
}

void G4UIcmdWith3Vector::SetParameterName
(const char * theNameX,const char * theNameY,const char * theNameZ,
G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamX = GetParameter(0);
  theParamX->SetParameterName(theNameX);
  theParamX->SetOmittable(omittable);
  theParamX->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamY = GetParameter(1);
  theParamY->SetParameterName(theNameY);
  theParamY->SetOmittable(omittable);
  theParamY->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZ = GetParameter(2);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
}

void G4UIcmdWith3Vector::SetDefaultValue(G4ThreeVector vec)
{
  G4UIparameter * theParamX = GetParameter(0);
  theParamX->SetDefaultValue(vec.x());
  G4UIparameter * theParamY = GetParameter(1);
  theParamY->SetDefaultValue(vec.y());
  G4UIparameter * theParamZ = GetParameter(2);
  theParamZ->SetDefaultValue(vec.z());
}

