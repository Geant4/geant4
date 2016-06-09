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
// $Id: G4PhysicalVolumeMassScene.cc,v 1.3 2004/11/11 16:06:49 johna Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

G4PhysicalVolumeMassScene::G4PhysicalVolumeMassScene ():
  fVolume (0.),
  fMass (0.),
  fpLastPV (0),
  fPVPCount (0),
  fLastDepth (0),
  fLastDensity (0.),
  fCurrentDepth (0),
  fpCurrentPV (0),
  fpCurrentLV (0)
{}

G4PhysicalVolumeMassScene::~G4PhysicalVolumeMassScene () {}

void G4PhysicalVolumeMassScene::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4PhysicalVolumeMassScene::Reset ()
{
  fVolume = 0.;
  fMass = 0.;
  fpLastPV = 0;
  fPVPCount = 0;
  fLastDepth = 0;
  fLastDensity = 0.;
  fDensityStack.clear();
  fCurrentDepth = 0;
  fpCurrentPV = 0;
  fpCurrentLV = 0;

}

void G4PhysicalVolumeMassScene::AccrueMass (const G4VSolid& solid)
{
  if (fpCurrentPV != fpLastPV) {
    fpLastPV = fpCurrentPV;
    fPVPCount = 0;
  }

  G4double currentVolume;
  G4double currentDensity;
  G4Polyhedron* pPolyhedron = solid.GetPolyhedron();
  if (pPolyhedron) {
    G4Material* pMaterial;
    G4VPVParameterisation* pP = fpCurrentPV->GetParameterisation();
    if (pP) {
      pMaterial = pP -> ComputeMaterial (fPVPCount++, fpCurrentPV);
    } else {
      pMaterial = fpCurrentLV->GetMaterial();
    }
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

  if (fCurrentDepth == 0) fVolume = currentVolume;

  if (fCurrentDepth > fLastDepth) {
    fDensityStack.push_back (fLastDensity);
  } else if (fCurrentDepth < fLastDepth) {
    fDensityStack.pop_back();
  }
  fLastDepth = fCurrentDepth;
  fLastDensity = currentDensity;
  G4double motherDensity = 0.;
  if (fCurrentDepth > 0) motherDensity = fDensityStack.back();

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
	   << fpCurrentPV->GetName() <<
      "\", copy "
	   << fpCurrentPV->GetCopyNo() <<
      ".  Larger than mother?"
	   << G4endl;
  }
}
