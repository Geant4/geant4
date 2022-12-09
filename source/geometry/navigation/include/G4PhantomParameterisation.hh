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
// class G4PhantomParameterisation
//
// Class description:
// 
// Describes regular parameterisations: a set of boxes of equal dimension
// in the x, y and z dimensions. The G4PVParameterised volume using this
// class must be placed inside a volume that is completely filled by these
// boxes.

// History:
// - Created: P.Arce, May 2007
//---------------------------------------------------------------------
#ifndef G4PhantomParameterisation_HH
#define G4PhantomParameterisation_HH

#include <vector>

#include "G4Types.hh"
#include "G4VPVParameterisation.hh"
#include "G4AffineTransform.hh"

class G4VPhysicalVolume;
class G4VTouchable; 
class G4VSolid;
class G4Material;

// Dummy forward declarations ...

class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Polycone;
class G4Polyhedra;

class G4PhantomParameterisation : public G4VPVParameterisation
{
  public:

    G4PhantomParameterisation();
   ~G4PhantomParameterisation();

    virtual void ComputeTransformation(const G4int, G4VPhysicalVolume *) const;
  
    virtual G4VSolid* ComputeSolid(const G4int, G4VPhysicalVolume *);
  
    virtual G4Material* ComputeMaterial(const G4int repNo, 
                                             G4VPhysicalVolume* currentVol,
                                       const G4VTouchable* parentTouch=nullptr);
  // Dummy declarations ...

    void ComputeDimensions (G4Box &, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&, const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&, const G4int,
                            const G4VPhysicalVolume*) const {}
  
    void BuildContainerSolid( G4VPhysicalVolume* pPhysicalVol );
    void BuildContainerSolid( G4VSolid* pMotherSolid );
      // Save as container solid the parent of the voxels. Check that the
      // voxels fill it completely.

    virtual G4int GetReplicaNo( const G4ThreeVector& localPoint,
                                const G4ThreeVector& localDir );
      // Get the voxel number corresponding to the point in the container
      // frame. Use 'localDir' to avoid precision problems at the surfaces.

    // Set and Get methods

    inline void SetMaterials(std::vector<G4Material*>& mates );

    inline void SetMaterialIndices( std::size_t* matInd );

    void SetVoxelDimensions( G4double halfx, G4double halfy, G4double halfz );
    void SetNoVoxels( std::size_t nx, std::size_t ny, std::size_t nz );
    
    inline G4double GetVoxelHalfX() const;
    inline G4double GetVoxelHalfY() const;
    inline G4double GetVoxelHalfZ() const;
    inline std::size_t GetNoVoxelsX() const;
    inline std::size_t GetNoVoxelsY() const;
    inline std::size_t GetNoVoxelsZ() const;
    inline std::size_t GetNoVoxels() const;

    inline std::vector<G4Material*> GetMaterials() const;
    inline std::size_t* GetMaterialIndices() const;
    inline G4VSolid* GetContainerSolid() const;

    G4ThreeVector GetTranslation(const G4int copyNo ) const;

    G4bool SkipEqualMaterials() const;
    void SetSkipEqualMaterials( G4bool skip );

    std::size_t GetMaterialIndex( std::size_t nx, std::size_t ny, std::size_t nz) const;
    std::size_t GetMaterialIndex( std::size_t copyNo) const;

    G4Material* GetMaterial( std::size_t nx, std::size_t ny, std::size_t nz) const;
    G4Material* GetMaterial( std::size_t copyNo ) const;

    void CheckVoxelsFillContainer( G4double contX, G4double contY,
                                   G4double contZ ) const;
      // Check that the voxels fill it completely.

  private:

    void ComputeVoxelIndices(const G4int copyNo, std::size_t& nx,
                                   std::size_t& ny, std::size_t& nz ) const;
      // Convert the copyNo to voxel numbers in x, y and z.

    void CheckCopyNo( const G4long copyNo ) const;
      // Check that the copy number is within limits.

  protected:

    G4double fVoxelHalfX = 0.0, fVoxelHalfY = 0.0, fVoxelHalfZ = 0.0;
      // Half dimension of voxels (assume they are boxes).
    std::size_t fNoVoxelsX = 0, fNoVoxelsY = 0, fNoVoxelsZ = 0;
      // Number of voxel in x, y and z dimensions.
    std::size_t fNoVoxelsXY = 0;
      // Number of voxels in x times number of voxels in y (for speed-up).
    std::size_t fNoVoxels = 0;
      // Total number of voxels (for speed-up).

    std::vector<G4Material*> fMaterials;
      // List of materials of the voxels.
    std::size_t* fMaterialIndices = nullptr;
      // Index in fMaterials that correspond to each voxel.

    G4VSolid* fContainerSolid = nullptr;
      // Save as container solid the parent of the voxels.
      // Check that the voxels fill it completely.

    G4double fContainerWallX=0.0, fContainerWallY=0.0, fContainerWallZ=0.0;
      // Save position of container wall for speed-up.

    G4double kCarTolerance;
      // Relative surface tolerance.

    G4bool bSkipEqualMaterials = true;
      // Flag to skip surface when two voxel have same material or not
};

#include "G4PhantomParameterisation.icc"

#endif
