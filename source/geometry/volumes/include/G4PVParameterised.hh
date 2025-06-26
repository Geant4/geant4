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
// G4PVParameterised
//
// Class description:
//
// Represents many touchable detector elements differing in their
// positioning and dimensions. Both are calculated by means
// of a G4VParameterisation object. The positioning is assumed to
// be dominant along a specified Cartesian axis.

// Author: Paul Kent (CERN), 29 July 1995 - first non-stub version
// ----------------------------------------------------------------------
#ifndef G4PVPARAMETERISED_HH
#define G4PVPARAMETERISED_HH 1

#include "G4PVReplica.hh"

/**
 * @brief G4PVParameterised represents many touchable detector elements
 * differing in their positioning and dimensions. Both are calculated by means
 * of a G4VParameterisation object. The positioning is assumed to be dominant
 * along a specified Cartesian axis.
 */

class G4PVParameterised : public G4PVReplica
{
  public:

    /**
     * Replicates the volume 'nReplicas' times using the paramaterisation
     * 'pParam', within the mother volume 'pMotherLogical'.
     * The positioning of the replicas is dominant along the specified axis.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the replica.
     *  @param[in] pMotherLogical Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the replication.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] pParam Pointer to the provided parameterisation algorithm.
     *  @param[in] pSurfChk Boolean flag, if true activates check for overlaps
     *             with existing volumes (false by default).
     */
    G4PVParameterised(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4LogicalVolume* pMotherLogical,
                      const EAxis pAxis,
                      const G4int nReplicas,
                            G4VPVParameterisation* pParam,
                            G4bool pSurfChk = false);
 
    /**
     * Similar to the constructor above, except for the mother pointer's type
     * being here a G4VPhysicalVolume.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the replica.
     *  @param[in] pMother Pointer to the physical volume of the mother.
     *  @param[in] pAxis The axis along which do the replication.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] pParam Pointer to the provided parameterisation algorithm.
     *  @param[in] pSurfChk Boolean flag, if true activates check for overlaps
     *             with existing volumes (false by default).
     */
    G4PVParameterised(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4VPhysicalVolume* pMother,
                      const EAxis pAxis,
                      const G4int nReplicas,
                            G4VPVParameterisation* pParam,
                            G4bool pSurfChk = false);

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4PVParameterised(__void__&);

    /**
     * Virtual empty Destructor.
     */
    ~G4PVParameterised() override;

    /**
     * Returns true to identify if it is a parameterised physical volume.
     */
    G4bool IsParameterised() const override;

    /**
     * Returns the volume type characterisation.
     */
    EVolume VolumeType() const final;
   
    /**
     * Returns the current pointer to the parameterisation algorithm.
     */
    G4VPVParameterisation* GetParameterisation() const override;

    /**
     * Fills arguments with the attributes from the base replica.
     * @param[in,out] axis Axis of parameterisation returned.
     * @param[in,out] nReplicas The number of replica copies.
     * @param[in,out] width Width of the replica object.
     * @param[in,out] offset Potential offset in replication.
     * @param[in,out] consuming Flag of replica characterisation (always false
     *                for parameterisations).
     */
    void GetReplicationData(EAxis& axis,
                            G4int& nReplicas,
                            G4double& width,
                            G4double& offset,
                            G4bool& consuming) const override;

    /**
     * Sets code and can prepare for special type of regular volumes.
     */
    void SetRegularStructureId( G4int code ) override; 

    /**
     * Verifies if each instance of the parameterised volume is overlapping
     * with other instances or with the mother volume. Provides default
     * resolution for the number of points to be generated and verified.
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
    G4bool CheckOverlaps(G4int res = 1000, G4double tol = 0.,
                         G4bool verbose = true, G4int maxErr = 1) override;

  private:

    /** The pointer to the parameterisation algorithm. */
    G4VPVParameterisation* fparam = nullptr; 
};

#endif
