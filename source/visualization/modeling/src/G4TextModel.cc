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
// $Id: G4TextModel.cc,v 1.4 2001-08-14 18:43:32 johna Exp $
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

void G4TextModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  sceneHandler.BeginPrimitives ();
  sceneHandler.AddPrimitive (fText);
  sceneHandler.EndPrimitives ();
}
