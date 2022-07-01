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
// John Allison  15th February 2019
// An artificial scene to reuse G4VScene code to calculate a bounding extent.

#include "G4BoundingExtentScene.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"

G4BoundingExtentScene::G4BoundingExtentScene (G4VModel* pModel)
:fpModel(pModel)
{}

G4BoundingExtentScene::~G4BoundingExtentScene () {}

void G4BoundingExtentScene::ProcessVolume(const G4VSolid& solid)
{
  G4VisExtent newExtent = solid.GetExtent ();
  if (fpCurrentObjectTransformation) {
    newExtent.Transform (*fpCurrentObjectTransformation);
  }
  AccrueBoundingExtent (newExtent);

  // Curtail descent - can assume daughters are contained within mother...
  G4PhysicalVolumeModel* pPVM = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVM) pPVM->CurtailDescent();
}

void G4BoundingExtentScene::ResetBoundingExtent ()
{
  fExtent = G4VisExtent::GetNullExtent();
  fpCurrentObjectTransformation = 0;
}

void G4BoundingExtentScene::AccrueBoundingExtent (const G4VisExtent& newExtent)
{

  if (fExtent == G4VisExtent::GetNullExtent()) {  // First time.

    fExtent = newExtent;

  } else {

    if (newExtent.GetXmin() < fExtent.GetXmin()) fExtent.SetXmin(newExtent.GetXmin());
    if (newExtent.GetYmin() < fExtent.GetYmin()) fExtent.SetYmin(newExtent.GetYmin());
    if (newExtent.GetZmin() < fExtent.GetZmin()) fExtent.SetZmin(newExtent.GetZmin());
    if (newExtent.GetXmax() > fExtent.GetXmax()) fExtent.SetXmax(newExtent.GetXmax());
    if (newExtent.GetYmax() > fExtent.GetYmax()) fExtent.SetYmax(newExtent.GetYmax());
    if (newExtent.GetZmax() > fExtent.GetZmax()) fExtent.SetZmax(newExtent.GetZmax());

  }
}
