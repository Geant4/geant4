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

#include "G4VtkPolydataPolylinePipeline.hh"

#include "G4VisAttributes.hh"

#include "vtkActor.h"
#include "vtkProperty.h"

std::size_t G4VtkPolydataPolylinePipeline::MakeHash(const G4VisAttributes* pVA)
{
  return std::hash<G4VisAttributes>{}(*pVA);
}

G4VtkPolydataPolylinePipeline::G4VtkPolydataPolylinePipeline(G4String nameIn,
                                                             const G4VtkVisContext& vcIn,
                                                             const G4VisAttributes* pVisAttributes)
  : G4VtkPolydataPipeline(nameIn, vcIn)
{
  // Set pipeline type
  SetTypeName(G4String("G4VtkPolydataPolylinePipeline"));

  // Get vis attributes
  G4Color colour = pVisAttributes->GetColour();
  G4double opacity = colour.GetAlpha();
  G4double lineWidth = pVisAttributes->GetLineWidth();

  // Setup actor and mapper
  GetActor()->GetProperty()->SetLineWidth(lineWidth);
  GetActor()->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
  GetActor()->GetProperty()->SetOpacity(opacity);
  GetActor()->SetVisibility(1);
}

void G4VtkPolydataPolylinePipeline::Print()
{
  G4cout << "G4VtkPolydataPolylinePipeline " << GetName() << G4endl;
  G4VtkPolydataPipeline::Print();
}
