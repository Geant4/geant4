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
// G4AssemblyVolume
//
// Class description:
//
// G4AssemblyVolume is a helper class to make the build process of geometry
// easier. It allows one to combine several volumes together in an arbitrary way
// in 3D space and then work with the result as with a single logical volume
// for placement.
// The resulting objects are independent copies of each of the assembled
// logical volumes. The placements are not, however, bound one to each other
// when placement is done. They are seen as independent physical volumes in
// space.

// Authors: R.Chytracek, J.Apostolakis, G.Cosmo (CERN), November 2000
// ----------------------------------------------------------------------
#ifndef G4_ASSEMBLYVOLUME_HH
#define G4_ASSEMBLYVOLUME_HH 1

#include <vector>

#include "G4Transform3D.hh"
#include "G4AssemblyTriplet.hh"

class G4VPhysicalVolume;

/**
 * @brief G4AssemblyVolume is a helper class to make the build process of
 * geometry easier. It allows one to combine several volumes together in an
 * arbitrary way in 3D space and then work with the result as with a single
 * logical volume for placement.
 * The resulting objects are independent copies of each of the assembled
 * logical volumes. The placements are not, however, bound one to each other
 * when placement is done. They are seen as independent physical volumes in
 * space.
 */

class G4AssemblyVolume
{
  public:

    /**
     * Default Constructor.
     */
    G4AssemblyVolume();    

    /**
     * Constructor.
     * The rotation matrix passed as argument can be nullptr (identity) or an
     * address even of an object on the upper stack frame. During assembly
     * imprint, a new matrix is created anyway and it is kept track of it so
     * it can be automatically deleted later at the end of the application.
     * This policy is adopted since user has no control on the way the
     * rotations are combined.
     *  @param[in] volume Pointer to the logical volume of the assembly.
     *  @param[in] translation Translation vector of the assembly.
     *  @param[in] rotation Pointer to the rotation matrix of the assembly.
     */
    G4AssemblyVolume( G4LogicalVolume* volume,
                      G4ThreeVector& translation,
                      G4RotationMatrix* rotation);
    
    /**
     * Destructor.
     * At destruction all the generated physical volumes and associated
     * rotation matrices of the imprints will be destroyed.
     */
    ~G4AssemblyVolume();

    /**
     * Places the given volume 'pPlacedVolume' inside the assembly.
     *
     * The adopted approach:
     *
     * - Place it w.r.t. the assembly coordinate system.
     *   This step is applied to each of the participating volumes.
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
     *
     * The rotation matrix passed as argument can be nullptr (identity) or an
     * address even of an object on the upper stack frame. During assembly
     * imprint, a new matrix is created anyway and it is kept track of it so
     * it can be automatically deleted later at the end of the application.
     * This policy is adopted since user has no control on the way the
     * rotations are combined.
     *  @param[in] pPlacedVolume Pointer to the logical volume of the volume
     *             to be added to the assembly.
     *  @param[in] translation Translation vector of the volume.
     *  @param[in] rotation Pointer to the rotation matrix of the volume.
     */
    void AddPlacedVolume( G4LogicalVolume* pPlacedVolume,
                          G4ThreeVector& translation,
                          G4RotationMatrix* rotation);

    /**
     * The same as previous method, except that it takes the complete 3D
     * transformation in space as its argument.
     *  @param[in] pPlacedVolume Pointer to the logical volume of the volume
     *             to be added to the assembly.
     *  @param[in] transformation The 3D transformation in space.
     */
    void AddPlacedVolume( G4LogicalVolume* pPlacedVolume,
                          G4Transform3D&   transformation);

    /**
     * The same as previous method, but takes an assembly volume as argument.
     *  @param[in] pAssembly Pointer to the assembly volume to be added
     *             to the assembly.
     *  @param[in] transformation The 3D transformation in space.
     */
    void AddPlacedAssembly( G4AssemblyVolume* pAssembly,
                            G4Transform3D&    transformation);

    /**
     * The same as previous method, but takes an assembly volume 
     * as its argument with translation and rotation.
     *  @param[in] pAssembly Pointer to the assembly volume to be added
     *             to the assembly.
     *  @param[in] translation Translation vector of the volume.
     *  @param[in] rotation Pointer to the rotation matrix of the volume.
     */
    void AddPlacedAssembly( G4AssemblyVolume* pAssembly,
                            G4ThreeVector& translation,
                            G4RotationMatrix* rotation);

    /**
     * Creates instance of an assembly volume inside the given mother volume.
     *  @param[in] pMotherLV Pointer to the logical volume of the assembly.
     *  @param[in] translationInMother Translation vector of the imprint.
     *  @param[in] pRotationInMother Pointer to the rotation of the imprint.
     *  @param[in] copyNumBase Optional index to assign to the imprint.
     *  @param[in] surfCheck Flag to enable overlaps checking while imprinting.
     */
    void MakeImprint( G4LogicalVolume* pMotherLV,
                      G4ThreeVector& translationInMother,
                      G4RotationMatrix* pRotationInMother,
                      G4int copyNumBase = 0,
                      G4bool surfCheck = false );

    /**
     * The same as previous Imprint() method, but takes a complete 3D
     * transformation in space as its argument.
     *  @param[in] pMotherLV Pointer to the logical volume of the assembly.
     *  @param[in] transformation The 3D transformation in space of the imprint.
     *  @param[in] copyNumBase Optional index to assign to the imprint.
     *  @param[in] surfCheck Flag to enable overlaps checking while imprinting.
     */
    void MakeImprint( G4LogicalVolume* pMotherLV,
                      G4Transform3D&   transformation,
                      G4int copyNumBase = 0,
                      G4bool surfCheck = false );

    /**
     * To access the physical volumes imprinted through an iterator.
     */
    inline std::vector<G4VPhysicalVolume*>::iterator GetVolumesIterator();

    /**
     * Returns the total number of imprinted volumes of the assembly.
     */
    inline std::size_t TotalImprintedVolumes() const;

    /**
     * To access the 3D transformation in space for each imprint, given the ID.
     */
    inline G4Transform3D& GetImprintTransformation(unsigned int imprintID);

    /**
     * To access the created triplets in the assembly through an iterator.
     */
    inline std::vector<G4AssemblyTriplet>::iterator GetTripletsIterator();

    /**
     * Returns the total number of triplets in the assembly.
     */
    inline std::size_t TotalTriplets() const;
  
    /**
     * Returns the number of currently made imprints.
     */
    inline unsigned int GetImprintsCount() const;

    /**
     * Returns the number of existing instances of G4AssemblyVolume class.
     */
    unsigned int GetInstanceCount() const;

    /**
     * Returns the instance number of the assembly.
     */
    inline unsigned int GetAssemblyID()    const;
  
  protected:

    inline void SetInstanceCount( unsigned int value );
    inline void SetAssemblyID( unsigned int value );
 
    void InstanceCountPlus();
    void InstanceCountMinus();

    // Internal counting mechanism, used to compute unique the names of
    // physical volumes created by MakeImprint() methods.
    //
    inline void SetImprintsCount( unsigned int value );
    inline void ImprintsCountPlus();
    inline void ImprintsCountMinus();

  private:    

    /**
     * Function for placement of the given assembly in the given mother
     * (called recursively if the assembly contains an assembly).
     */
    void MakeImprint( G4AssemblyVolume* pAssembly,
                      G4LogicalVolume*  pMotherLV,
                      G4Transform3D&    transformation,
                      G4int copyNumBase = 0,
                      G4bool surfCheck = false );

  private:

    /** Participating volumes represented as a vector of
     *  <logical volume, translation, rotation>. */
    std::vector<G4AssemblyTriplet> fTriplets;

    /** We need to keep list of physical volumes created by MakeImprint()
     * in order to be able to cleanup the objects when not needed anymore.
     * This requires the user to keep assembly objects in memory during the
     * whole job or during the life-time of G4Navigator, logical volume store
     * and physical volume store keep pointers to physical volumes generated
     * by the assembly volume.
     * When an assembly object is about to die it will destroy all its
     * generated physical volumes and rotation matrices as well ! */
    std::vector<G4VPhysicalVolume*> fPVStore;

    /** Number of imprints of the given assembly volume. */
    unsigned int fImprintsCounter;

    /** Class instance counter. */
    static G4ThreadLocal unsigned int fsInstanceCounter;

    /** Assembly object ID derived from instance counter at construction time. */
    unsigned int fAssemblyID = 0;

    /** Container of transformations for each imprint (used in GDML). */
    std::map<unsigned int, G4Transform3D> fImprintsTransf;
};

#include "G4AssemblyVolume.icc"

#endif
