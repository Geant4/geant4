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
//
// 
// John Allison  3rd April 2001
// Model which knows how to draw text.

#include "G4TextModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"

#include "G4UnitsTable.hh"
#include <sstream>

G4TextModel::~G4TextModel () {}

G4TextModel::G4TextModel (const G4Text& g4Text, const G4Transform3D& transform)
: fG4Text(g4Text)
{
  fType = "G4TextModel";
  std::ostringstream oss;
  oss << "G4TextModel: \"" << fG4Text.GetText()
      << "\" at " << G4BestUnit(g4Text.GetPosition(),"Length")
      << "with size " << g4Text.GetScreenSize()
      << " with offsets " << g4Text.GetXOffset() << ',' << g4Text.GetYOffset();
  fGlobalTag = oss.str();
  fGlobalDescription = fGlobalTag;

  fG4Text.SetPosition(fG4Text.GetPosition().transform(transform));
}

void G4TextModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  sceneHandler.BeginPrimitives ();
  sceneHandler.AddPrimitive (fG4Text);
  sceneHandler.EndPrimitives ();
}
