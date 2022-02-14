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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
// This class developed by John Allison, March 2021
//

#include "ICRP110PhantomPseudoScene.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

// Given a world volume, fill maps provided by the instantiator
ICRP110PhantomPseudoScene::ICRP110PhantomPseudoScene
(G4PhysicalVolumeModel* pvm  // input...the following are outputs
 ,std::multimap<const G4Material*,G4ThreeVector>& positionByMaterial
 ,std::map     <const G4Material*,G4Colour>&      colourByMaterial)
: fpPVModel(pvm)
, fPositionByMaterial(positionByMaterial)
, fColourByMaterial  (colourByMaterial)
{}

ICRP110PhantomPseudoScene::~ICRP110PhantomPseudoScene()
{}

void ICRP110PhantomPseudoScene::AddSolid(const G4Box&)
{
  const G4Material* material = fpPVModel->GetCurrentLV()->GetMaterial();
  const G4Colour& colour = fpPVModel->GetCurrentLV()->GetVisAttributes()->GetColour();
  const G4ThreeVector& position = fpCurrentObjectTransformation->getTranslation();

  // Fill maps for vis action
  fPositionByMaterial.insert(std::make_pair(material,position));
  fColourByMaterial[material] = colour;
}

void ICRP110PhantomPseudoScene::ProcessVolume(const G4VSolid&)
{
  G4ExceptionDescription ed;
  ed << "ICRP110PhantomPseudoScene::ProcessVolume called for solid that is NOT a G4Box.";
  G4Exception("ICRP110PhantomPseudoScene::ProcessVolume", "ICRP110Phantom_0000",
	      FatalException, ed, "Contact @allison");
}
