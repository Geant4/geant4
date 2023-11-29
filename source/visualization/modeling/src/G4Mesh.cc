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
// G4Mesh captures and validates a parameterisation, which we
// call a "mesh". This is typically intended for meshes with
// a large number of parameterisations, such as a medical phantom.
//
// G4Mesh is used by G4PhysicalVolumeModel if and only if
// G4ModelingParameters::fSpecialMeshRendering is set and if the
// name matches one in G4ModelingParameters::fSpecialMeshVolumes,
// if any. Then, if a valid mesh is found it calls the overriding
// implementation of G4VGraphicsScene::AddCompound(const G4Mesh&).
//
// To set the above parameters use the following commands in the
// standard Geant4 Visualisation System:
//    /vis/viewer/set/specialMeshRendering
//    /vis/viewer/set/specialMeshRenderingOption
//    /vis/viewer/set/specialMeshVolumes
// See guidance on the above commmands for more detail.
//
// Note that if no special mesh volumes are specified,
// G4PhysicalVolumeModel will test all volumes, and therefore
// it will capture *all* parameterisations. This is not usually
// a problem, since there is usually only one, but to be
// selective you have to /vis/viewer/set/specialMeshVolumes.
//
// The specified G4VPhysicalVolume is searched for a
// parameterisation. If none is found it will have a type "invalid"
// and it should simply be destroyed (as in G4PhysicalVolumeModel).
// The overhead of an invalid attempt is small.

#include "G4Mesh.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4VNestedParameterisation.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Tet.hh"

std::map<G4int,G4String> G4Mesh::fEnumMap = {
  {invalid,"invalid"},
  {rectangle,"rectangle"},
  {nested3DRectangular,"nested3Drectangular"},
  {cylinder,"cylinder"},
  {sphere,"sphere"},
  {tetrahedron,"tetrahedron"}
};

G4Mesh::G4Mesh (G4VPhysicalVolume* containerVolume,const G4Transform3D& transform)
: fpContainerVolume(containerVolume)
, fpParameterisedVolume(nullptr)
, fMeshType(invalid)
, fMeshDepth(0)
, fTransform(transform)
{
  if (fpContainerVolume == nullptr) return;
    
  G4VPhysicalVolume* pv0 = fpContainerVolume;
  G4VPhysicalVolume* pv1 = nullptr;
  G4VPhysicalVolume* pv2 = nullptr;
  G4VPhysicalVolume* pv3 = nullptr;
  G4LogicalVolume*   lv0 = pv0->GetLogicalVolume();
  G4LogicalVolume*   lv1 = nullptr;
  G4LogicalVolume*   lv2 = nullptr;

  // Check if this is a container for a parameterisation.
  // A simple parameterisation may only be one level.
  // Nested parameterisations may be 2- or 3-level.
  G4bool isContainer = false;
  if (lv0->GetNoDaughters()) {
    fMeshDepth++;
    pv1 = lv0->GetDaughter(0);
    lv1 = pv1->GetLogicalVolume();
    if (dynamic_cast<G4PVParameterised*>(pv1)) {
      isContainer = true;
      fpParameterisedVolume = pv1;
    } else if (lv1->GetNoDaughters()) {
      fMeshDepth++;
      pv2 = lv1->GetDaughter(0);
      lv2 = pv2->GetLogicalVolume();
      if (dynamic_cast<G4PVParameterised*>(pv2) &&
          dynamic_cast<G4VNestedParameterisation*>(pv2->GetParameterisation())) {
        isContainer = true;
        fpParameterisedVolume = pv2;
      } else if (lv2->GetNoDaughters()) {
        fMeshDepth++;
        pv3 = lv2->GetDaughter(0);
        if (dynamic_cast<G4PVParameterised*>(pv3) &&
            dynamic_cast<G4VNestedParameterisation*>(pv3->GetParameterisation())) {
          isContainer = true;
          fpParameterisedVolume = pv3;
        }
      }
    }
  }

  if (isContainer) {

    // Get type
    G4VSolid* pEndSol = fpParameterisedVolume->GetLogicalVolume()->GetSolid ();
    if (dynamic_cast<G4Box*>(pEndSol)) {
      fMeshType = rectangle;
      auto pBox = static_cast<G4Box*>(pEndSol);
      f3DRPs.fHalfX = pBox->GetXHalfLength();
      f3DRPs.fHalfY = pBox->GetYHalfLength();
      f3DRPs.fHalfZ = pBox->GetZHalfLength();
    } else if (dynamic_cast<G4Tet*>(pEndSol)) {
      fMeshType = tetrahedron;
    } else if (dynamic_cast<G4Tubs*>(pEndSol)) {
      fMeshType = cylinder;
    } else if (dynamic_cast<G4Sphere*>(pEndSol)) {
      fMeshType = sphere;
    }
    
    // Special case for rectangular nested paramaterisation - extra information
    if (fMeshDepth == 3 && fMeshType == rectangle) {
      auto nestedParam3 = dynamic_cast<G4VNestedParameterisation*>(pv3);
      if (nestedParam3) {
        fMeshType = nested3DRectangular;
        pv1->GetReplicationData
        (f3DRPs.fAxis1,f3DRPs.fNreplica1,f3DRPs.fWidth1,f3DRPs.fOffset1,f3DRPs.fConsuming1);
        pv2->GetReplicationData
        (f3DRPs.fAxis2,f3DRPs.fNreplica2,f3DRPs.fWidth2,f3DRPs.fOffset2,f3DRPs.fConsuming2);
        pv3->GetReplicationData
        (f3DRPs.fAxis3,f3DRPs.fNreplica3,f3DRPs.fWidth3,f3DRPs.fOffset3,f3DRPs.fConsuming3);
      }
    }
  }
}

G4Mesh::~G4Mesh () {}

std::ostream& operator << (std::ostream& os, const G4Mesh& mesh) {
  os << "G4Mesh: ";
  os << "\nContainer: " << mesh.GetContainerVolume()->GetName();
  const auto& map = mesh.GetEnumMap();
  const auto& typeEntry = map.find(mesh.GetMeshType());
  G4String type;
  if (typeEntry != map.end()) {
    type = typeEntry->second;
  } else {
    type = "unrecognised";
  }
  os << "\nType: " << type;
  os << "\nDepth: " << mesh.GetMeshDepth();
  os << "\nTranslation: " << mesh.GetTransform().getTranslation();
  os << "\nRotation: " << mesh.GetTransform().getRotation();
  if (mesh.GetMeshType() == G4Mesh::rectangle &&
      mesh.GetMeshDepth() == 3) {
    // Print ThreeDRectangleParameters
  }
  return os;
}
