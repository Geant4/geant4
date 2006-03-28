//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PhysicalVolumeMassScene.cc,v 1.5 2006-03-28 16:46:27 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

void G4PhysicalVolumeMassScene::AccrueMass (const G4VSolid& solid)
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
  G4double currentDensity = pCurrentMaterial->GetDensity();
  /* Using G4Polyhedron... (gives slightly different answers on Tubs, e.g.).
  G4Polyhedron* pPolyhedron = solid.GetPolyhedron();
  if (pPolyhedron) {
    G4Material* pMaterial;
    G4VPVParameterisation* pP = pCurrentPV->GetParameterisation();
    if (pP) {
      pMaterial = pP -> ComputeMaterial (fPVPCount++, pCurrentPV);
    } else {
      pMaterial = pCurrentLV->GetMaterial();
    }
    assert(pMaterial == pCurrentMaterial);
    currentVolume = pPolyhedron->GetVolume();
    currentDensity = pMaterial->GetDensity();
  } else {
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
    G4cout <<
      "G4PhysicalVolumeMassScene::AccrueMass: WARNING:"
      "\n  Mass going negative for \""
	   << pCurrentPV->GetName() <<
      "\", copy "
	   << pCurrentPV->GetCopyNo() <<
      ".  Larger than mother?"
	   << G4endl;
  }
}
