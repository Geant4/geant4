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
// Class Description:
// G4Mesh encapsulates and validates a nested parameterisation, which we
// call a "mesh". If a valid mesh cannot be created out of this
// G4VPhysicalVolume* (which will probably be most common), it will
// have a type "invalid". Then, usually, it may simply be destroyed.
// The overhead of an invalid attempt is expected to be small.
// Class Description - End:

#ifndef G4MESH_HH
#define G4MESH_HH

#include "G4Transform3D.hh"

class G4VPhysicalVolume;

class G4Mesh
{
 public:  // With description
  enum MeshType
  {
    invalid,
    rectangle,
    cylinder,
    sphere
  };

  G4Mesh(G4VPhysicalVolume* containerVolume, const G4Transform3D&);
  virtual ~G4Mesh();

  const std::map<G4int,G4String>& GetEnumMap() const { return fEnumMap; }
  G4VPhysicalVolume* GetContainerVolume() const { return fpContainerVolume; }
  MeshType GetMeshType() const { return fMeshType; }
  G4int GetMeshDepth() const { return fMeshDepth; }
  const G4Transform3D& GetTransform() const { return fTransform; }

 private:
  static std::map<G4int,G4String> fEnumMap;
  G4VPhysicalVolume* fpContainerVolume;
  MeshType fMeshType;
  G4int fMeshDepth;
  G4Transform3D fTransform;
};

std::ostream& operator << (std::ostream& os, const G4Mesh& mesh);

#endif
