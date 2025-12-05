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
// G4PVDivision
//
// Class description:
//
// Represents many touchable detector elements differing only in their
// positioning. The elements' positions are calculated by means of a simple
// linear formula.

// Author: Pedro Arce (CIEMAT), 09.05.2001 - Initial version
// ----------------------------------------------------------------------
#ifndef G4PVDIVISION_HH
#define G4PVDIVISION_HH 1

#include "geomdefs.hh"
#include "G4PVReplica.hh"
#include "G4VDivisionParameterisation.hh"

class G4LogicalVolume;
class G4VSolid;

/**
 * @brief G4PVDivision represents many touchable detector elements differing
 * only in their positioning. The elements' positions are calculated by means
 * of a simple linear formula.
 * 
 * Division may occur along:
 *
 * o Cartesian axes (kXAxis,kYAxis,kZAxis)
 *
 *   The divisions, of specified width have coordinates of
 *   form (-width*(nReplicas-1)*0.5+n*width,0,0) where n=0.. nReplicas-1
 *   for the case of kXAxis, and are unrotated.
 *
 * o Radial axis (cylindrical polar) (kRho)
 *
 *   The divisions are cons/tubs sections, centred on the origin
 *   and are unrotated.
 *   They have radii of width*n+offset to width*(n+1)+offset
 *                      where n=0..nReplicas-1
 *
 * o Phi axis (cylindrical polar) (kPhi)
 *   The divisions are `phi sections' or wedges, and of cons/tubs form
 *   They have phi of offset+n*width to offset+(n+1)*width where
 *   n=0..nReplicas-1
 */

class G4PVDivision : public G4PVReplica
{
  public:
    
    /**
     * Constructor with number of divisions and width.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] width The witdh of the divided slice along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     */
    G4PVDivision(const G4String& pName,
                       G4LogicalVolume* pLogical,
                       G4LogicalVolume* pMother,
                 const EAxis pAxis,
                 const G4int nReplicas,
                 const G4double width,
                 const G4double offset );

    /**
     * Constructor with number of divisions.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] offset The optional offset distance from mother's border.
     */
    G4PVDivision(const G4String& pName,
                       G4LogicalVolume* pLogical,
                       G4LogicalVolume* pMother,
                 const EAxis pAxis,
                 const G4int nReplicas,
                 const G4double offset );

    /**
     * Constructor with width.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] width The witdh of the divided slice along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     */
    G4PVDivision(const G4String& pName,
                       G4LogicalVolume* pLogical,
                       G4LogicalVolume* pMother,
                 const EAxis pAxis,
                 const G4double width,
                 const G4double offset );

    /**
     * Constructor in mother physical volume (same as first constructor).
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the division.
     *  @param[in] pMother Pointer to the physical volume of the mother.
     *  @param[in] pAxis The axis along which do the division.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] width The witdh of the divided slice along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     */
    G4PVDivision(const G4String& pName,
                       G4LogicalVolume* pLogical,
                       G4VPhysicalVolume* pMother,
                 const EAxis pAxis,
                 const G4int nReplicas,
                 const G4double width,
                 const G4double offset);

    /**
     * Destructor.
     */
    ~G4PVDivision() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4PVDivision(const G4PVDivision&) = delete;
    G4PVDivision& operator=(const G4PVDivision&) = delete;

    /**
     * Not used.
     */
    G4bool IsMany() const override;

    /**
     * Returns true.
     */
    G4bool IsReplicated() const override;

    /**
     * Returns the number of divisions.
     */
    G4int GetMultiplicity() const override;

    /**
     * Returns true to identify if it is a parameterised physical volume.
     */
    G4bool IsParameterised() const override;

    /**
     * Returns the pointer to the parameterisation algorithm.
     */
    G4VPVParameterisation* GetParameterisation() const override;

    /**
     * Fills arguments with the attributes from the base replica.
     * @param[in,out] axis Axis of parameterisation returned.
     * @param[in,out] nReplicas The number of division copies.
     * @param[in,out] width Width of the division slice.
     * @param[in,out] offset Potential offset in replication.
     * @param[in,out] consuming Flag of replica characterisation (always false
     *                for parameterisations).
     */
    void GetReplicationData( EAxis& axis,
                             G4int& nReplicas,
                             G4double& width,
                             G4double& offset,
                             G4bool& consuming ) const override;

    /**
     * Returns the axis along which the division is made.
     *  @returns The string Id of the axis.
     */
    EAxis GetDivisionAxis() const;
   
    /**
     * Returns the volume type characterisation.
     *  @returns The string Id of the volume type, i.e. 'kParameterised'.
     */
    EVolume VolumeType() const override;

    /**
     * Methods to identify volumes that can have revised 'regular' navigation.
     * Currently divisions do not qualify for this.
     */
    G4bool IsRegularStructure() const override;    // Returns false
    G4int GetRegularStructureId() const override;  // Returns 0

  private:

    /**
     * Invoked in constructors to check and set the parameters for the
     * parameterisation.
     */
    void CheckAndSetParameters( const EAxis pAxis,
                                const G4int nDivs,
                                const G4double width,
                                const G4double offset, 
                                      DivisionType divType,
                                const G4LogicalVolume* pMotherLogical );

    /**
     * Sets the parameterisation according to the allowed type of solid
     * to be divided.
     */
    void SetParameterisation( G4LogicalVolume* motherLogical,
                        const EAxis pAxis,
                        const G4int nReplicas,
                        const G4double width,
                        const G4double offset, 
                              DivisionType divType );

    /**
     * Logging the axis when an error in division occurs.
     */
    void ErrorInAxis( EAxis axis, G4VSolid* solid );

  protected:

    /** Axis of optimisation. */
    EAxis faxis;

    /** Axis of division. */
    EAxis fdivAxis;

    /** Number of divisions. */
    G4int fnReplicas = 0;

    /** Width of the division slice and potential offset. */
    G4double fwidth = 0.0, foffset = 0.0;

    /** Pointer to the parameterisation algorithm. */
    G4VDivisionParameterisation* fparam = nullptr; 
};

#endif
