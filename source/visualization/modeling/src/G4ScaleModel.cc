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
// $Id: G4ScaleModel.cc,v 1.2 2001-07-24 21:50:09 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  21st July 2001.
// Model which knows how to draw scale.

#include "G4ScaleModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"

G4ScaleModel::~G4ScaleModel () {}

G4ScaleModel::G4ScaleModel (const G4Scale& scale): fScale(scale) {
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
  sceneHandler.AddPrimitive (fScale);
}
