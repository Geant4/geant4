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
// G4VExternalPhysicalVolume
//
// Class description:
//
// Base class to represent a physical volume managed by an external 
// sub-navigator.
// 
// Intial assumptions:
//   * volume type is similar to G4PVPlacement -- not replicated
//   * external navigator may provide 'many'/Boolean operation

// Author: John Apostolakis (CERN), October 2019
// ----------------------------------------------------------------------
#ifndef G4VEXTERNALPHYSICSVOLUME_HH
#define G4VEXTERNALPHYSICSVOLUME_HH 1

#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"

/**
 * @brief G4VExternalPhysicalVolume is a base class to represent a physical
 * volume managed by an external sub-navigator.
 */

class G4VExternalPhysicalVolume : public G4VPhysicalVolume
{
  public:

    /**
     * Constructor implementing G4VPhysicalVolume signature.
     *  @param[in] pRot The pointer to the rotation in the mother volume.
     *  @param[in] tlate Traslation vector in the mother volume.
     *  @param[in] pCurrentLogical Pointer to its logical volume.
     *  @param[in] pName The volume name.
     *  @param[in] pMotherLogical Pointer to the logical volume of the mother.
     */
    G4VExternalPhysicalVolume( G4RotationMatrix* pRot,
                               const G4ThreeVector& tlate,
                               G4LogicalVolume* pCurrentLogical,                            
                               const G4String& pName,
                               G4VPhysicalVolume* pMother );
       
    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4VExternalPhysicalVolume(__void__&);

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4VExternalPhysicalVolume(const G4VExternalPhysicalVolume&) = delete;
    G4VExternalPhysicalVolume& operator=(const G4VExternalPhysicalVolume&) = delete;

    /**
     * Default Destructor.
     */
    ~G4VExternalPhysicalVolume() override;

    /**
     * Verifies if the placed volume is overlapping with existing daughters
     * or with the mother volume. Provides default resolution for the number
     * of points to be generated and verified.
     * A tolerance for the precision of the overlap check can be specified,
     * by default it is set to maximum precision.
     * Reports a maximum of overlaps errors according to parameter in input.
     *  @param[in] res The number of points to generate on volume's surface.
     *  @param[in] tol The precision tolerance for the overlap check, below
     *             which to ignore overlaps (deafult is maximim precision).
     *  @param[in] verbose Verbosity mode (default is true).
     *  @param[in] maxErr Maximum of overlaps errors to report (default is 1).
     *  @returns True if an overlap occurs.
     */
    G4bool CheckOverlaps(G4int res=1000, G4double tol=0.,
                         G4bool verbose=true, G4int maxErr=1) override = 0;

    /**
     * Returns the volume type characterisation, i.e. kExternal.
     */
    EVolume VolumeType() const final;

    /**
     * Stub methods, unused for placed volumes.
     */
    G4bool IsMany() const final;
    void SetMany(G4bool overlap);
    G4bool IsReplicated() const final;
    G4bool IsParameterised() const final;
    G4VPVParameterisation* GetParameterisation() const final;
    void GetReplicationData(EAxis& axis,
                            G4int& nReplicas,
                            G4double& width,
                            G4double& offset,
                            G4bool& consuming) const final;

    /**
     * Used only for specialised repeated volumes. Always false with Id 0.
     */
    G4bool  IsRegularStructure() const final; 
    G4int   GetRegularStructureId() const final; 

  private:

    G4bool fMany = false;
};

#endif
