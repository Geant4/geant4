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
// $Id: G4PhantomParameterisation.hh 73435 2013-08-27 11:06:54Z gcosmo $
//
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
// - Created.    P. Arce, May 2007
// *********************************************************************

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
  public:  // with description

    G4PhantomParameterisation();
   ~G4PhantomParameterisation();

    virtual void ComputeTransformation(const G4int, G4VPhysicalVolume *) const;
  
    virtual G4VSolid* ComputeSolid(const G4int, G4VPhysicalVolume *);
  
    virtual G4Material* ComputeMaterial(const G4int repNo, 
                                              G4VPhysicalVolume *currentVol,
                                        const G4VTouchable *parentTouch=0);
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
  
    void BuildContainerSolid( G4VPhysicalVolume *pPhysicalVol );
    void BuildContainerSolid( G4VSolid *pMotherSolid );
      // Save as container solid the parent of the voxels. Check that the
      // voxels fill it completely.

    virtual G4int GetReplicaNo( const G4ThreeVector& localPoint,
                        const G4ThreeVector& localDir );
      // Get the voxel number corresponding to the point in the container
      // frame. Use 'localDir' to avoid precision problems at the surfaces.

    // Set and Get methods

    inline void SetMaterials(std::vector<G4Material*>& mates );

    inline void SetMaterialIndices( size_t* matInd );

    void SetVoxelDimensions( G4double halfx, G4double halfy, G4double halfz );
    void SetNoVoxel( size_t nx, size_t ny, size_t nz );
    
    inline G4double GetVoxelHalfX() const;
    inline G4double GetVoxelHalfY() const;
    inline G4double GetVoxelHalfZ() const;
    inline size_t GetNoVoxelX() const;
    inline size_t GetNoVoxelY() const;
    inline size_t GetNoVoxelZ() const;
    inline size_t GetNoVoxel() const;

    inline std::vector<G4Material*> GetMaterials() const;
    inline size_t* GetMaterialIndices() const;
    inline G4VSolid* GetContainerSolid() const;

    G4ThreeVector GetTranslation(const G4int copyNo ) const;

    G4bool SkipEqualMaterials() const;
    void SetSkipEqualMaterials( G4bool skip );

    size_t GetMaterialIndex( size_t nx, size_t ny, size_t nz) const;
    size_t GetMaterialIndex( size_t copyNo) const;

    G4Material* GetMaterial( size_t nx, size_t ny, size_t nz) const;
    G4Material* GetMaterial( size_t copyNo ) const;

    void CheckVoxelsFillContainer( G4double contX, G4double contY,
                                   G4double contZ ) const;
      // Check that the voxels fill it completely.

  private:

    void ComputeVoxelIndices(const G4int copyNo, size_t& nx,
                                   size_t& ny, size_t& nz ) const;
      // Convert the copyNo to voxel numbers in x, y and z.

    void CheckCopyNo( const G4int copyNo ) const;
      // Check that the copy number is within limits.

  protected:

    G4double fVoxelHalfX,fVoxelHalfY,fVoxelHalfZ;
      // Half dimension of voxels (assume they are boxes).
    size_t fNoVoxelX,fNoVoxelY,fNoVoxelZ;
      // Number of voxel in x, y and z dimensions.
    size_t fNoVoxelXY;
      // Number of voxels in x times number of voxels in y (for speed-up).
    size_t fNoVoxel;
      // Total number of voxels (for speed-up).

    std::vector<G4Material*> fMaterials;
      // List of materials of the voxels.
    size_t* fMaterialIndices;
      // Index in fMaterials that correspond to each voxel.

    G4VSolid* fContainerSolid;
      // Save as container solid the parent of the voxels.
      // Check that the voxels fill it completely.

    G4double fContainerWallX, fContainerWallY, fContainerWallZ;
      // Save position of container wall for speed-up.

    G4double kCarTolerance;
      // Relative surface tolerance.

    G4bool bSkipEqualMaterials;
      // Flag to skip surface when two voxel have same material or not
};

#include "G4PhantomParameterisation.icc"

#endif
