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

#include "G4VtkTextPipeline.hh"

#include "G4Text.hh"
#include "G4VisAttributes.hh"
#include "G4VtkViewer.hh"

#include "vtkBillboardTextActor3D.h"
#include "vtkTextProperty.h"

std::size_t G4VtkTextPipeline::MakeHash(const G4Text& text, const G4VtkVisContext& /*vc*/,
                                        const G4VisAttributes* pVA)
{
  std::size_t hash = std::hash<G4VisAttributes>{}(*pVA);

  G4double x = text.GetPosition().x();
  G4double y = text.GetPosition().y();
  G4double z = text.GetPosition().z();

  std::hash_combine(hash, std::hash<G4String>{}(text.GetText()));
  std::hash_combine(hash, std::hash<G4double>{}(x));
  std::hash_combine(hash, std::hash<G4double>{}(y));
  std::hash_combine(hash, std::hash<G4double>{}(z));

  return hash;
}

G4VtkTextPipeline::G4VtkTextPipeline(const G4Text& text, const G4VtkVisContext& vc,
                                     const G4VisAttributes* pVA)
  : G4VVtkPipeline(text.GetText().c_str(), "G4VtkTextPipeline", vc, false, vc.fViewer->renderer)
{
  G4double x = text.GetPosition().x();
  G4double y = text.GetPosition().y();
  G4double z = text.GetPosition().z();

  G4Color colour = pVA->GetColour();
  G4double opacity = colour.GetAlpha();

  actor = vtkSmartPointer<vtkBillboardTextActor3D>::New();
  actor->SetInput(text.GetText().c_str());
  actor->SetPosition(x, y, z);
  actor->GetTextProperty()->SetFontSize(text.GetScreenSize());
  actor->GetTextProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
  actor->GetTextProperty()->SetOpacity(opacity);

  vc.fViewer->renderer->AddActor(actor);
}

void G4VtkTextPipeline::Enable()
{
  actor->SetVisibility(1);
}

void G4VtkTextPipeline::Disable()
{
  actor->SetVisibility(0);
}

void G4VtkTextPipeline::Print()
{
  G4VVtkPipeline::Print();
}

void G4VtkTextPipeline::Modified()
{
  G4VVtkPipeline::Modified();
}

void G4VtkTextPipeline::Clear()
{
  renderer->RemoveActor(actor);
  G4VVtkPipeline::Clear();
}

void G4VtkTextPipeline::SetText(const G4String& text)
{
  actor->SetInput(text.c_str());
}
