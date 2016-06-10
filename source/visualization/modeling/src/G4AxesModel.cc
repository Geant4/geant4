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
// $Id: G4AxesModel.cc 83403 2014-08-21 15:07:30Z gcosmo $
//
// 
// John Allison  3rd April 2001
// Model which knows how to draw axes.

#include "G4AxesModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4VisAttributes.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4ArrowModel.hh"
#include "G4TextModel.hh"

G4AxesModel::~G4AxesModel ()
{
  delete fZAnnotationModel;
  delete fZLabelModel;
  delete fZAxisModel;
  delete fYAnnotationModel;
  delete fYLabelModel;
  delete fYAxisModel;
  delete fXAnnotationModel;
  delete fXLabelModel;
  delete fXAxisModel;
}

G4AxesModel::G4AxesModel
(G4double x0, G4double y0, G4double z0, G4double length,
 G4double arrowWidth, const G4String& colourString,
 const G4String& description,
 G4bool withAnnotation,
 G4double textSize):
  fXAxisModel(0),
  fXLabelModel(0),
  fXAnnotationModel(0),
  fYAxisModel(0),
  fYLabelModel(0),
  fYAnnotationModel(0),
  fZAxisModel(0),
  fZLabelModel(0),
  fZAnnotationModel(0)
{
  fType = "G4AxesModel";
  fGlobalTag = fType;
  fGlobalDescription = fType + ": " + description;
  fExtent = G4VisExtent
    (x0, x0+length, y0, y0+length, z0, z0+length);

  G4Colour colour(1,1,1,1);  // Default white and opaque (unless "auto").
  G4bool autoColour = false;
  if (colourString == "auto") autoColour = true;
  else {
    if (!G4Colour::GetColour(colourString, colour)) {
      G4ExceptionDescription ed;
      ed << "Colour \"" << colourString
	 << "\" not found.  Defaulting to white and opaque.";
      G4Exception
	("G4AxesModel::G4AxesModel",
	 "modeling0012", JustWarning, ed);
    }
  }

  G4String annotation = G4BestUnit(length,"Length");

  G4Text* text = 0;
  G4VisAttributes* va = 0;

  G4Colour xColour(colour);
  if (autoColour) xColour = G4Colour::Red();
  fXAxisModel = new G4ArrowModel
    (x0, y0, z0, x0+length, y0, z0, arrowWidth,
     xColour, "x-axis: " + description);
  if (withAnnotation) {
    text = new G4Text("x",G4Point3D(x0+1.05*length, y0, z0));
    text->SetScreenSize(textSize);
    text->SetOffset(0.5*textSize,0.5*textSize);
    text->SetLayout(G4Text::centre);
    va = new G4VisAttributes(xColour);
    text->SetVisAttributes(va);
    fXLabelModel = new G4TextModel(*text);
    delete text;
    text = new G4Text(annotation,G4Point3D(x0+0.8*length, y0, z0));
    text->SetScreenSize(textSize);
    text->SetOffset(-1.5*textSize,-1.5*textSize);
    text->SetLayout(G4Text::centre);
    va = new G4VisAttributes(xColour);
    text->SetVisAttributes(va);
    fXAnnotationModel = new G4TextModel(*text);
    delete text;
  }

  G4Colour yColour(colour);
  if (autoColour) yColour = G4Colour::Green();
  fYAxisModel = new G4ArrowModel
    (x0, y0, z0, x0, y0+length, z0, arrowWidth,
     yColour, "y-axis: " + description);
  if (withAnnotation) {
    text = new G4Text("y",G4Point3D(x0, y0+1.05*length, z0));
    text->SetScreenSize(textSize);
    text->SetOffset(0.5*textSize,0.5*textSize);
    text->SetLayout(G4Text::centre);
    va = new G4VisAttributes(yColour);
    text->SetVisAttributes(va);
    fYLabelModel = new G4TextModel(*text);
    delete text;
    text = new G4Text(annotation,G4Point3D(x0, y0+0.8*length, z0));
    text->SetScreenSize(textSize);
    text->SetOffset(-1.5*textSize,-1.5*textSize);
    text->SetLayout(G4Text::centre);
    va = new G4VisAttributes(yColour);
    text->SetVisAttributes(va);
    fYAnnotationModel = new G4TextModel(*text);
    delete text;
  }

  G4Colour zColour(colour);
  if (autoColour) zColour = G4Colour::Blue();
  fZAxisModel = new G4ArrowModel
    (x0, y0, z0, x0, y0, z0+length, arrowWidth,
     zColour, "z-axis: " + description);
  if (withAnnotation) {
    text = new G4Text("z",G4Point3D(x0, y0, z0+1.05*length));
    text->SetScreenSize(textSize);
    text->SetOffset(0.5*textSize,0.5*textSize);
    text->SetLayout(G4Text::centre);
    va = new G4VisAttributes(zColour);
    text->SetVisAttributes(va);
    fZLabelModel = new G4TextModel(*text);
    delete text;
    text = new G4Text(annotation,G4Point3D(x0, y0, z0+0.8*length));
    text->SetScreenSize(textSize);
    text->SetOffset(-1.5*textSize,-1.5*textSize);
    text->SetLayout(G4Text::centre);
    va = new G4VisAttributes(zColour);
    text->SetVisAttributes(va);
    fZAnnotationModel = new G4TextModel(*text);
    delete text;
  }
}

void G4AxesModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  if (fXAxisModel)       fXAxisModel->DescribeYourselfTo(sceneHandler);
  if (fXLabelModel)      fXLabelModel->DescribeYourselfTo(sceneHandler);
  if (fXAnnotationModel) fXAnnotationModel->DescribeYourselfTo(sceneHandler);
  if (fYAxisModel)       fYAxisModel->DescribeYourselfTo(sceneHandler);
  if (fYLabelModel)      fYLabelModel->DescribeYourselfTo(sceneHandler);
  if (fYAnnotationModel) fYAnnotationModel->DescribeYourselfTo(sceneHandler);
  if (fZAxisModel)       fZAxisModel->DescribeYourselfTo(sceneHandler);
  if (fZLabelModel)      fZLabelModel->DescribeYourselfTo(sceneHandler);
  if (fZAnnotationModel) fZAnnotationModel->DescribeYourselfTo(sceneHandler);
}
