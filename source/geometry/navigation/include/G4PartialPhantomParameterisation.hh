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
// $Id: G4PartialPhantomParameterisation.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4PartialPhantomParameterisation
//
// Class description:
// 
// Describes partial regular parameterisations, i.e. the voxels do not 
// completely fill the container in the three dimensions

// History:
// - Created.    P. Arce, September 2010
// *********************************************************************

#ifndef G4PartialPhantomParameterisation_HH
#define G4PartialPhantomParameterisation_HH

#include <vector>
#include <set>
#include <map>

#include "G4Types.hh"
#include "G4PhantomParameterisation.hh"
#include "G4AffineTransform.hh"

class G4VPhysicalVolume;
class G4VTouchable; 
class G4VSolid;
class G4Material;

class G4PartialPhantomParameterisation : public G4PhantomParameterisation
{
  public:  // with description

    G4PartialPhantomParameterisation();
   ~G4PartialPhantomParameterisation();

    void ComputeTransformation(const G4int, G4VPhysicalVolume *) const;
  
    G4Material* ComputeMaterial(const G4int repNo, 
                                      G4VPhysicalVolume *currentVol,
                                const G4VTouchable *parentTouch=0);

    G4int GetReplicaNo( const G4ThreeVector& localPoint,
                        const G4ThreeVector& localDir );
      // Get the voxel number corresponding to the point in the container
      // frame. Use 'localDir' to avoid precision problems at the surfaces.

    G4ThreeVector GetTranslation(const G4int copyNo ) const;

    size_t GetMaterialIndex( size_t nx, size_t ny, size_t nz) const;
    size_t GetMaterialIndex( size_t copyNo) const;

    G4Material* GetMaterial( size_t nx, size_t ny, size_t nz) const;
    G4Material* GetMaterial( size_t copyNo ) const;

    void SetFilledIDs(   std::multimap<G4int,G4int> fid ){
      fFilledIDs = fid; 
    }

    void SetFilledMins( std::map< G4int, std::map<G4int,G4int> > fmins ) {
      fFilledMins = fmins;
    }

    void BuildContainerWalls();

  private:

    void ComputeVoxelIndices(const G4int copyNo, size_t& nx,
                                   size_t& ny, size_t& nz ) const;
      // Convert the copyNo to voxel numbers in x, y and z.

    void CheckCopyNo( const G4int copyNo ) const;
      // Check that the copy number is within limits.

  private:

    std::multimap<G4int,G4int> fFilledIDs;
    std::map< G4int, std::map<G4int,G4int> > fFilledMins;
};

#endif
