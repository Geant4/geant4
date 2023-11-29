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

#include "G4VtkPolydataSpherePipeline.hh"

#include "G4VisAttributes.hh"
#include "G4VtkViewer.hh"
#include "G4VtkVisContext.hh"

#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkVertexGlyphFilter.h"

std::size_t G4VtkPolydataSpherePipeline::MakeHash(const G4VisAttributes* pVA)
{
  return std::hash<G4VisAttributes>{}(*pVA);
}

G4VtkPolydataSpherePipeline::G4VtkPolydataSpherePipeline(G4String nameIn, const G4VtkVisContext& vc,
                                                         const G4VisAttributes* pVisAttributes)
  : G4VtkPolydataPipeline(nameIn, vc)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataSpherePipeline"));

  // Get vis attributes
  G4Color colour = pVisAttributes->GetColour();
  G4double opacity = colour.GetAlpha();
  G4double lineWidth = pVisAttributes->GetLineWidth();

  // vertex glyph filter
  auto glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputData(polydata);
  AddFilter(glyphFilter);

  // Setup actor and mapper
  mapper->SetInputConnection(glyphFilter->GetOutputPort());

  GetActor()->GetProperty()->SetLineWidth(lineWidth);
  GetActor()->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
  GetActor()->GetProperty()->SetOpacity(opacity);
  GetActor()->SetVisibility(1);

  GetActor()->GetProperty()->SetRenderPointsAsSpheres(true);
  GetActor()->GetProperty()->SetPointSize(vc.fSize * 5);

  vc.fViewer->renderer->AddActor(GetActor());
}

void G4VtkPolydataSpherePipeline::Print()
{
  G4cout << "G4VtkPolydataSpherePipeline " << GetName() << G4endl;
  G4VtkPolydataPipeline::Print();
}
