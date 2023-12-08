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

#include "G4VtkPolydataPolyline2DPipeline.hh"

#include "G4Polyline.hh"
#include "G4VisAttributes.hh"
#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkCoordinate.h"
#include "vtkLine.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"

std::size_t G4VtkPolydataPolyline2DPipeline::MakeHash(const G4VisAttributes* pVA)
{
  return std::hash<G4VisAttributes>{}(*pVA);
}

G4VtkPolydataPolyline2DPipeline::G4VtkPolydataPolyline2DPipeline(
  G4String nameIn, const G4VtkVisContext& vcIn, const G4VisAttributes* pVisAttributes)
  : G4VtkPolydataPipeline(nameIn, vcIn)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataPolyline2DPipeline"));

  // Get vis attributes
  G4Color colour = pVisAttributes->GetColour();
  G4double opacity = colour.GetAlpha();
  G4double lineWidth = pVisAttributes->GetLineWidth();

  // make 3d actor invisible and remove from renderer
  GetActor()->SetVisibility(0);
  vc.fViewer->renderer->RemoveActor(GetActor());

  // Setup 2D mapper
  mapper2D = vtkSmartPointer<vtkPolyDataMapper2D>::New();
  mapper2D->SetInputConnection(GetFinalFilter()->GetOutputPort());
  vtkSmartPointer<vtkCoordinate> coordinate = vtkSmartPointer<vtkCoordinate>::New();
  coordinate->SetCoordinateSystemToNormalizedViewport();
  mapper2D->SetTransformCoordinate(coordinate);

  // Setup 2D actor
  actor2D = vtkSmartPointer<vtkActor2D>::New();
  GetActor2D()->GetProperty()->SetLineWidth(lineWidth);
  GetActor2D()->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
  GetActor2D()->GetProperty()->SetOpacity(opacity);
  GetActor2D()->SetMapper(mapper2D);
  GetActor2D()->SetVisibility(1);
  // GetActor2D()->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();

  // add to renderer
  vc.fViewer->renderer->AddActor(GetActor2D());
}

void G4VtkPolydataPolyline2DPipeline::Clear()
{
  renderer->RemoveActor(GetActor2D());
  G4VtkPolydataPipeline::Clear();
}

void G4VtkPolydataPolyline2DPipeline::Modified()
{
  G4VtkPolydataPipeline::Modified();
  GetActor2D()->Modified();
}

void G4VtkPolydataPolyline2DPipeline::SetPolydata(const G4Polyline& polyline)
{
  // Data data
  const size_t nLines = polyline.size();

  for (size_t i = 0; i < nLines; ++i) {
    auto id =
      polydataPoints->InsertNextPoint((polyline[i].x() + 1) / 2.0, (polyline[i].y() + 1) / 2.0, 0);

    if (i < nLines - 1) {
      vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
      line->GetPointIds()->SetId(0, id);
      line->GetPointIds()->SetId(1, id + 1);
      polydataCells->InsertNextCell(line);
    }
  }
}
