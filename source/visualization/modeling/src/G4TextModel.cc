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
// $Id: G4TextModel.cc,v 1.6 2005/06/07 16:53:40 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

G4TextModel::G4TextModel (const G4Text& g4Text): fG4Text(g4Text) {
  std::ostringstream oss;
  oss << "G4TextModel: \"" << fG4Text.GetText()
      << "\" at " << G4BestUnit(g4Text.GetPosition(),"Length")
      << "with size " << g4Text.GetScreenSize()
      << " with offsets " << g4Text.GetXOffset() << ',' << g4Text.GetYOffset();
  fGlobalTag = oss.str();
  fGlobalDescription = fGlobalTag;
}

void G4TextModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  sceneHandler.BeginPrimitives ();
  sceneHandler.AddPrimitive (fG4Text);
  sceneHandler.EndPrimitives ();
}
