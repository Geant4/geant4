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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
// Code based on the extended example DICOM
//
#ifndef ICRP110PhantomNestedParameterisation_HH
#define ICRP110PhantomNestedParameterisation_HH

#include <vector>
#include <map>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VTouchable.hh"
#include "G4VNestedParameterisation.hh" 

class G4VPhysicalVolume;
class G4VSolid;
class G4Material;
class G4VisAttributes;

// CSG Entities which may be parameterised/replicated
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

class ICRP110PhantomNestedParameterisation : public G4VNestedParameterisation
{
public:
  
  explicit ICRP110PhantomNestedParameterisation(const G4ThreeVector& voxelSize,
                                     std::vector<G4Material*>& mat,
                                     G4int fnX_ = 0, G4int fnY_ = 0, G4int fnZ_ = 0);
                                     // the total number of voxels along X, Y and Z 
                                     // are initialised to zero 
                                     
  ~ICRP110PhantomNestedParameterisation(); 
  
  
  virtual G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
                              const G4int repNo, 
                              const G4VTouchable *parentTouch );
  G4int       GetNumberOfMaterials() const;
  G4Material* GetMaterial(G4int idx) const;

  G4int GetMaterialIndex( G4int copyNo) const;
  
  void SetMaterialIndices( size_t* matInd ){ fMaterialIndices = matInd;}
  // This method passes the information of the matID associated to each voxel
  // from the DetectorConstruction to the NestedParameterisation class
  
  void SetNoVoxel( G4int nx, G4int ny, G4int nz ); 
  // This method passes the total number of voxels along X, Y and Z from 
  // the DetectorConstruction to the NestedParameterisation class
  
  void ComputeTransformation(const G4int no,
                             G4VPhysicalVolume *currentPV) const;
  
  // Additional standard Parameterisation methods, 
  // which can be optionally defined, in case solid is used.
  void ComputeDimensions(G4Box &, const G4int,
                         const G4VPhysicalVolume *) const;
  
  
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
  
  void ReadColourData();
  
  using G4VNestedParameterisation::ComputeMaterial;
  
private:
  
  G4double fdX,fdY,fdZ; // Half of the voxels along X, Y and Z 
  G4int fnX,fnY,fnZ;    // Number of voxels along X, Y and Z
  std::vector<G4Material*> fMaterials; // Vector with materials
  size_t* fMaterialIndices; // Index of the material associated to each voxel
  std::map<G4String,G4VisAttributes*> fColours;
};
#endif
