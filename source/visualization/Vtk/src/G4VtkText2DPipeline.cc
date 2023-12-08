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
// Created by Stewart Boogert on 05/03/2023.
//

#include "G4VtkText2DPipeline.hh"

#include "G4Text.hh"
#include "G4VisAttributes.hh"
#include "G4VtkViewer.hh"

#include "vtkTextActor.h"
#include "vtkTextProperty.h"

std::size_t G4VtkText2DPipeline::MakeHash(const G4Text& text, const G4VtkVisContext& /*vc*/,
                                          const G4VisAttributes* pVA)
{
  std::size_t hash = std::hash<G4VisAttributes>{}(*pVA);

  G4double x = text.GetPosition().x();
  G4double y = text.GetPosition().y();

  std::hash_combine(hash, std::hash<G4String>{}(text.GetText()));
  std::hash_combine(hash, std::hash<G4double>{}(x));
  std::hash_combine(hash, std::hash<G4double>{}(y));

  return hash;
}

G4VtkText2DPipeline::G4VtkText2DPipeline(const G4Text& text, const G4VtkVisContext& vcIn,
                                         const G4VisAttributes* pVA)
  : G4VVtkPipeline(text.GetText().c_str(), "G4VtkText2DPipeline", vcIn, false, vcIn.fViewer->renderer)
{
  G4double x = text.GetPosition().x();
  G4double y = text.GetPosition().y();

  G4Color colour = pVA->GetColour();
  G4double opacity = colour.GetAlpha();

  actor = vtkSmartPointer<vtkTextActor>::New();
  actor->SetInput(text.GetText().c_str());
  actor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  // actor->SetTextScaleModeToViewport();
  actor->SetPosition((x + 1.) / 2.0, (y + 1.) / 2.);
  actor->GetTextProperty()->SetFontSize(text.GetScreenSize());
  actor->GetTextProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
  actor->GetTextProperty()->SetOpacity(opacity);

  switch (text.GetLayout()) {
    case G4Text::Layout::left:
      actor->GetTextProperty()->SetJustificationToLeft();
      break;
    case G4Text::Layout::centre:
      actor->GetTextProperty()->SetJustificationToCentered();
      break;
    case G4Text::Layout::right:
      actor->GetTextProperty()->SetJustificationToRight();
      break;
  }

  vc.fViewer->renderer->AddActor(actor);
}

void G4VtkText2DPipeline::Enable()
{
  actor->SetVisibility(1);
}

void G4VtkText2DPipeline::Disable()
{
  actor->SetVisibility(0);
}

void G4VtkText2DPipeline::Print()
{
  G4VVtkPipeline::Print();
}

void G4VtkText2DPipeline::Modified()
{
  G4VVtkPipeline::Modified();
}

void G4VtkText2DPipeline::Clear()
{
  renderer->RemoveViewProp(actor);

  G4VVtkPipeline::Clear();
}

void G4VtkText2DPipeline::SetText(const G4String& text)
{
  actor->SetInput(text.c_str());
}
