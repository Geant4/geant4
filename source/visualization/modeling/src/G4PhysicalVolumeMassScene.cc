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
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#include "G4PhysicalVolumeMassScene.hh"

#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedron.hh"
#include "G4Material.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"

#define G4warn G4cout

G4PhysicalVolumeMassScene::G4PhysicalVolumeMassScene
(G4PhysicalVolumeModel* pPVModel):
  fpPVModel (pPVModel),
  fVolume (0.),
  fMass (0.),
  fpLastPV (0),
  fPVPCount (0),
  fLastDepth (0),
  fLastDensity (0.)
{}

G4PhysicalVolumeMassScene::~G4PhysicalVolumeMassScene () {}

void G4PhysicalVolumeMassScene::Reset ()
{
  fVolume = 0.;
  fMass = 0.;
  fpLastPV = 0;
  fPVPCount = 0;
  fLastDepth = 0;
  fLastDensity = 0.;
  fDensityStack.clear();
}

void G4PhysicalVolumeMassScene::ProcessVolume (const G4VSolid& solid)
{
  G4int currentDepth = fpPVModel->GetCurrentDepth();
  G4VPhysicalVolume* pCurrentPV = fpPVModel->GetCurrentPV();
  //G4LogicalVolume* pCurrentLV = fpPVModel->GetCurrentLV();
  G4Material* pCurrentMaterial = fpPVModel->GetCurrentMaterial();

  if (pCurrentPV != fpLastPV) {
    fpLastPV = pCurrentPV;
    fPVPCount = 0;
  }

  G4double currentVolume = ((G4VSolid&)solid).GetCubicVolume();
  G4double currentDensity = pCurrentMaterial? pCurrentMaterial->GetDensity() : 0.;
  /* Using G4Polyhedron... (gives slightly different answers on Tubs, e.g.).
  G4Polyhedron* pPolyhedron = solid.GetPolyhedron();
  if (!pPolyhedron) {
    G4cout << 
      "G4PhysicalVolumeMassScene::AccrueMass: WARNING:"
      "\n  No G4Polyhedron for" << solid.GetEntityType() <<
      ".  \"" << solid.GetName() << "\" will not be accounted."
      "\n  It will be as though not there, i.e., the density as its mother."
      "\n  Its daughters will still be found and accounted."
	   << G4endl;
    currentVolume = 0.;
    currentDensity = 0.;
  }
  */

  if (currentDepth == 0) fVolume = currentVolume;

  if (currentDepth > fLastDepth) {
    fDensityStack.push_back (fLastDensity);
  } else if (currentDepth < fLastDepth) {
    fDensityStack.pop_back();
  }
  fLastDepth = currentDepth;
  fLastDensity = currentDensity;
  G4double motherDensity = 0.;
  if (currentDepth > 0) motherDensity = fDensityStack.back();

  G4double subtractedMass = currentVolume * motherDensity;
  G4double addedMass = currentVolume * currentDensity;
  fMass -= subtractedMass;
  fMass += addedMass;
  /* Debug
  G4cout << "current vol = "
	 << G4BestUnit (currentVolume,"Volume")
	 << ", current density = "
	 << G4BestUnit (currentDensity, "Volumic Mass")
	 << ", mother density = "
	 << G4BestUnit (motherDensity, "Volumic Mass")
	 << G4endl;
  G4cout << "Subtracted mass = " << G4BestUnit (subtractedMass, "Mass")
	 << ", added mass = " << G4BestUnit (addedMass, "Mass")
	 << ", new mass = " << G4BestUnit (fMass, "Mass")
	 << G4endl;
  */
  if (fMass < 0.) {
    G4warn <<
      "G4PhysicalVolumeMassScene::AccrueMass: WARNING:"
      "\n  Mass going negative for \""
	   << pCurrentPV->GetName() <<
      "\", copy "
	   << pCurrentPV->GetCopyNo() <<
      ".  Larger than mother?"
	   << G4endl;
  }
}
