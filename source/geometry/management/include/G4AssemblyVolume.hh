#ifndef G4ASSEMBLYVOLUME_H
#define G4ASSEMBLYVOLUME_H 1

#include "G4AssemblyTriplet.hh"
#include "g4std/vector"

/**
 * Class:       G4AssemblyVolume
 * Description: This is a helper class to make the build process of geometry
 *              easier. It allows to combine several volumes together in an
 *              arbitrary way in 3D space and then work with the result as with
 *              single logical volume for placement.
 *              The resulting objects are independent copies of each of the
 *              assembled logical volumes. The placements are not, however,
 *              bound one to each other when placement is done.
 *              They are seen as independent physical volumes in space.
 *              This is the main difference when comapred to the boolean solids
 *              or parametrized physical volumes.
 * 
 * Author:      Radovan Chytracek, John Apostolakis, Gabriele Cosmo
 * Version:     1.0
 * Date:        November 2000
 */
class G4AssemblyVolume {
public:

  G4AssemblyVolume();    
  G4AssemblyVolume( G4LogicalVolume* volume
                   ,G4ThreeVector& translation
                   ,G4RotationMatrix* rotation);
  ~G4AssemblyVolume();

  /**
   * Place the given volume inside the assembly.
   * 
   * The taken approach:
   *
   * - Place it w.r.t. the virtual coordinate system.
   *   This step is aplied to each of the participating volumes.
   *
   * The other possible approaches:
   *
   * - Place w.r.t. the firstly added volume.
   *   When placed the first, the virtual coordinate system becomes
   *   the coordinate system of the first one.
   *   Every next volume being added into the assembly will be placed
   *   w.r.t to the first one.
   *
   * - Place w.r.t the last placed volume.
   *   When placed the first, the virtual coordinate system becomes
   *   the coordinate system of the first one.
   *   Every next volume being added into the assembly will be placed
   *   w.r.t to the previous one.
   */
  void AddPlacedVolume( G4LogicalVolume* pPlacedVolume
                       ,G4ThreeVector& translation
                       ,G4RotationMatrix* rotation
                      );

  /**
   * Creates instance of an assembly volume inside the given mother volume
   */
  void MakeImprint( G4LogicalVolume* pMotherLV
                   ,G4ThreeVector& translationInMother
                   ,G4RotationMatrix* pRotationInMother
                  );

private:    

  /**
   * Participating volumes represented as a vector triplets of
   * < logical volume, rotation, translation >
   */
  G4std::vector<G4AssemblyTriplet> fTriplets; 
 
};

#endif //G4ASSEMBLYVOLUME_H

