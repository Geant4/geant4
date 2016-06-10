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
// $Id: G4ScaleModel.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison  21st July 2001.
// Model which knows how to draw scale.

#include "G4ScaleModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"

G4ScaleModel::~G4ScaleModel () {}

G4ScaleModel::G4ScaleModel (const G4Scale& scale): fScale(scale) {
  fType = "G4ScaleModel";
  fGlobalTag = "G4ScaleModel: " + fScale.GetAnnotation();
  switch (fScale.GetDirection()) {
  case G4Scale::x:
    fGlobalTag += " x";
    break;
  case G4Scale::y:
    fGlobalTag += " y";
    break;
  case G4Scale::z:
    fGlobalTag += " z";
    break;
  }
  fGlobalDescription = fGlobalTag;
}

void G4ScaleModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  sceneHandler.BeginPrimitives ();
  sceneHandler.AddPrimitive (fScale);
  sceneHandler.EndPrimitives ();
}
