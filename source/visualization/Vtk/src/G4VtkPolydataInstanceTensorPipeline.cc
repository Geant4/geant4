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

#include "G4VtkPolydataInstanceTensorPipeline.hh"

#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkTensorGlyphColor.h"

std::size_t G4VtkPolydataInstanceTensorPipeline::MakeHash(const G4Polyhedron& polyhedron,
                                                          const G4VtkVisContext& vc)
{
  // Get view parameters that the user can force through the vis attributes, thereby over-riding the
  // current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = vc.fDrawingStyle;

  std::size_t vhash = 0;
  std::hash_combine(vhash, static_cast<int>(drawing_style));

  std::size_t phash = std::hash<G4Polyhedron>{}(polyhedron);

  std::size_t hash = 0;
  std::hash_combine(hash, phash);
  std::hash_combine(hash, vhash);

  return hash;
}

G4VtkPolydataInstanceTensorPipeline::G4VtkPolydataInstanceTensorPipeline(G4String nameIn,
                                                                         const G4VtkVisContext& vc)
  : G4VtkPolydataInstancePipeline(nameIn, vc)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataInstanceTensorPipeline"));

  // make polydata from instance data
  instancePolydata = vtkSmartPointer<vtkPolyData>::New();
  instancePolydata->SetPoints(instancePosition);
  instancePolydata->GetPointData()->SetTensors(instanceTransform);
  instancePolydata->GetPointData()->SetVectors(instanceColour);
  instancePolydata->GetPointData()->SetScalars(instanceColour);

  // tensor glyph
  instanceTensorGlyph = vtkSmartPointer<vtkTensorGlyphColor>::New();
  instanceTensorGlyph->SetInputData(instancePolydata);
  instanceTensorGlyph->SetSourceConnection(GetFinalFilter()->GetOutputPort());
  instanceTensorGlyph->ColorGlyphsOn();
  instanceTensorGlyph->ScalingOff();
  instanceTensorGlyph->ThreeGlyphsOff();
  instanceTensorGlyph->ExtractEigenvaluesOff();
  instanceTensorGlyph->SetColorModeToScalars();
  AddFilter(instanceTensorGlyph);

  // set polydata mapper
  mapper->SetInputConnection(GetFinalFilter()->GetOutputPort());

  // set actor
  actor->SetMapper(mapper);
  actor->SetVisibility(1);

  // shading parameters
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.1);
  actor->GetProperty()->SetSpecularPower(1);

  // add to renderer
  vc.fViewer->renderer->AddActor(GetActor());
}

void G4VtkPolydataInstanceTensorPipeline::Print()
{
  G4cout << "G4VtkPolydataInstanceTensorPipeline " << GetName() << G4endl;
  G4VtkPolydataPipeline::Print();
}