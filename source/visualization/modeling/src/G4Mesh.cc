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
// John Allison  May 2021
//
// G4Mesh encapsulates and validates a nested parameterisation, which we
// call a "mesh". If a valid mesh cannot be created out of this
// G4VPhysicalVolume* (which will probably be most common), it will
// have a type "invalid". Then, usually, it may simply be destroyed.
// The overhead of an invalid attempt is expected to be small.

#include "G4Mesh.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4VNestedParameterisation.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"

std::map<G4int,G4String> G4Mesh::fEnumMap;

G4Mesh::G4Mesh (G4VPhysicalVolume* containerVolume,const G4Transform3D& transform)
: fpContainerVolume(containerVolume)
, fMeshType(invalid)
, fMeshDepth(0)
, fTransform(transform)
{
  if (fpContainerVolume == nullptr) return;

  static G4bool first = true;
  if (first) {
    first = false;
    fEnumMap[invalid]   = "invalid";
    fEnumMap[rectangle] = "rectangle";
    fEnumMap[cylinder]  = "cylinder";
    fEnumMap[sphere]    = "sphere";
  }

  const G4LogicalVolume* pLV = fpContainerVolume->GetLogicalVolume();

  // check if this is a container for a nested parameterisation
  G4bool isContainer = false;
  if (pLV->GetNoDaughters()) {
    fMeshDepth++;
    const auto d0 = pLV->GetDaughter(0);
    const auto d0LV = d0->GetLogicalVolume();
    if (d0LV->GetNoDaughters()) {
      fMeshDepth++;
      const auto d00 = d0LV->GetDaughter(0);
      const auto pvParam00 = dynamic_cast<G4PVParameterised*>(d00);
      if (pvParam00) {
	const auto param00 = pvParam00->GetParameterisation();
	const auto nestedParam00 = dynamic_cast<G4VNestedParameterisation*>(param00);
	if (nestedParam00) { // 2-deep mesh
	  isContainer = true;
	}
      } else {
	const auto d00LV = d00->GetLogicalVolume();
	if (d00LV->GetNoDaughters()) {
	  fMeshDepth++;
	  const auto d000 = d00LV->GetDaughter(0);
	  const auto pvParam000 = dynamic_cast<G4PVParameterised*>(d000);
	  if (pvParam000) {
	    const auto param000 = pvParam000->GetParameterisation();
	    const auto nestedParam000 = dynamic_cast<G4VNestedParameterisation*>(param000);
	    if (nestedParam000) {  // 3-deep mesh
	      isContainer = true;
	    }
	  }
	}
      }
    }
  }
  if (isContainer) {
    // Get type
    G4VSolid* pSol = pLV -> GetSolid ();
    if (dynamic_cast<G4Box*>(pSol)) {
      fMeshType = rectangle;
    } else if (dynamic_cast<G4Tubs*>(pSol)) {
      fMeshType = cylinder;
    } else if (dynamic_cast<G4Sphere*>(pSol)) {
      fMeshType = sphere;
    }
  }
}

G4Mesh::~G4Mesh () {}

std::ostream& operator << (std::ostream& os, const G4Mesh& mesh) {
  os << "G4Mesh: ";
  os << "\nContainer: " << mesh.GetContainerVolume()->GetName();
  os << "\nType: " << mesh.GetEnumMap().find(mesh.GetMeshType())->second;
  os << "\nDepth: " << mesh.GetMeshDepth();
  os << "\nTranslation: " << mesh.GetTransform().getTranslation();
  os << "\nRotation: " << mesh.GetTransform().getRotation();
  return os;
}
