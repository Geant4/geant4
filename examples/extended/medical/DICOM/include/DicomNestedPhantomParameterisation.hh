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
/// \file medical/DICOM/include/DicomNestedPhantomParameterisation.hh
/// \brief Definition of the DicomNestedPhantomParameterisation class
//

#ifndef DICOMNESTEDPARAMETERISATION_HH
#define DICOMNESTEDPARAMETERISATION_HH

#include <vector>
#include <map>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VNestedParameterisation.hh"

class G4VPhysicalVolume;
class G4VTouchable;
class G4VSolid;
class G4Material;
class G4VisAttributes;

// CSG Entities which may be parameterised/replicated
//
class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Ellipsoid;
class G4Orb;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

/// Implements a G4VNestedParameterisation

class DicomNestedPhantomParameterisation : public G4VNestedParameterisation
{
public:
    typedef std::map<G4String, G4VisAttributes*> ColourMap_t;

    static G4String defaultColorFile;

public:

  DicomNestedPhantomParameterisation(const G4ThreeVector& voxelSize,
                std::vector<G4Material*>& mat,
                G4int fnZ_ = 0, G4int fnY_ = 0, G4int fnX_ = 0,
                G4String colorFile = defaultColorFile);
  ~DicomNestedPhantomParameterisation();

  G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
                              const G4int repNo,
                              const G4VTouchable *parentTouch );
  // Must cope with parentTouch for navigator's SetupHierarchy

  G4int       GetNumberOfMaterials() const;
  G4Material* GetMaterial(G4int idx) const;
  // Needed to define materials for instances of Nested Parameterisation
  // Current convention: each call should return the materials
  // of all instances with the same mother/ancestor volume

  //unsigned int GetMaterialIndex( unsigned int nx, unsigned int ny,
  //                               unsigned int nz) const;
  unsigned int GetMaterialIndex( unsigned int copyNo) const;
  void SetMaterialIndices( size_t* matInd ) { fMaterialIndices = matInd; }
  void SetNoVoxel( unsigned int nx, unsigned int ny, unsigned int nz );

  void ComputeTransformation(const G4int no,
                             G4VPhysicalVolume *currentPV) const;

  // Additional standard Parameterisation methods,
  // which can be optionally defined, in case solid is used.
  void ComputeDimensions(G4Box &, const G4int,
                         const G4VPhysicalVolume *) const;

  const ColourMap_t& GetColourMap() const { return fColours; }
  ColourMap_t& GetColourMap() { return fColours; }

private:  // Dummy declarations to get rid of warnings ...

  void ComputeDimensions (G4Trd&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Trap&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Cons&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Sphere&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Ellipsoid&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Orb&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Torus&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Para&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Hype&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Tubs&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polycone&, const G4int,
                          const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polyhedra&, const G4int,
                          const G4VPhysicalVolume*) const {}

  void ReadColourData(G4String);

  using G4VNestedParameterisation::ComputeMaterial;

private:

  G4double fdX,fdY,fdZ;
  G4int fnX,fnY,fnZ;
  std::vector<G4Material*> fMaterials;
  size_t* fMaterialIndices; // Index in materials corresponding to each voxel
  ColourMap_t fColours;
  std::map<G4int, G4VisAttributes*> mColours;
  std::vector<G4double> fpZ;
};

#endif
