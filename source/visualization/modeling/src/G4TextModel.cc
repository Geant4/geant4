// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TextModel.cc,v 1.1 2001-04-11 13:39:34 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  3rd April 2001
// Model which knows how to draw text.

#include "G4TextModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"

G4TextModel::~G4TextModel () {}

G4TextModel::G4TextModel (const G4Text& text): fText(text) {
  fGlobalTag = "G4TextModel: " + fText.GetText();
  fGlobalDescription = fGlobalTag;
}

void G4TextModel::DescribeYourselfTo (G4VGraphicsScene& scene) {
  scene.AddPrimitive (fText);
}
