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
// G4PVReplica
//
// Class description:
//
// Represents many touchable detector elements differing only in their
// positioning. The elements' positions are calculated by means of a simple
// linear formula, and the elements completely fill the containing mother
// volume.

// Original author: Paul Kent (CERN), 29.07.1995 - First non-stub version
// - G.Cosmo, A.Dotti (CERN), 13.01.2013 - Modified for MT thread-safety
// ----------------------------------------------------------------------
#ifndef G4PVREPLICA_HH
#define G4PVREPLICA_HH

#include "G4VPhysicalVolume.hh"
#include "G4GeomSplitter.hh"

class G4ReplicaData
{
  /**
   * @brief G4ReplicaData encapsulates the fields of the class G4PVReplica that
   * may not be read-only. G4PVReplica inherits from the class G4VPhysicalVolume.
   * The fields from the ancestor that may not be read-only are handled
   * by the ancestor class.
   */

  public:

    void initialize() {}

    G4int fcopyNo = -1;
};

/** Implementation detail for use of G4ReplicaData objects. */
using G4PVRManager = G4GeomSplitter<G4ReplicaData>;

/**
 * @brief G4PVReplica represents many touchable detector elements differing
 * only in their positioning. The elements' positions are calculated by means
 * of a simple linear formula, and the elements completely fill the containing
 * mother volume.
 * 
 * Replication may occur along:
 *
 * o Cartesian axes (kXAxis,kYAxis,kZAxis)
 *
 *   The replications, of specified width have coordinates of
 *   form (-width*(nReplicas-1)*0.5+n*width,0,0) where n=0.. nReplicas-1
 *   for the case of kXAxis, and are unrotated.
 *
 * o Radial axis (cylindrical polar) (kRho)
 *
 *   The replications are cons/tubs sections, centred on the origin
 *   and are unrotated.
 *   They have radii of width*n+offset to width*(n+1)+offset
 *                      where n=0..nReplicas-1
 *
 * o Phi axis (cylindrical polar) (kPhi)
 *   The replications are `phi sections' or wedges, and of cons/tubs form
 *   They have phi of offset+n*width to offset+(n+1)*width where
 *   n=0..nReplicas-1
 */

class G4PVReplica : public G4VPhysicalVolume
{
  public:

    /**
     * Replicates the volume 'nReplicas' times along the specified axis
     * within the mother volume 'pMother' and filling completely the mother.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the replica.
     *  @param[in] pMother Pointer to the logical volume of the mother.
     *  @param[in] pAxis The axis along which do the replication.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] width The witdh of the replicated object along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     */
    G4PVReplica(const G4String& pName,
                      G4LogicalVolume* pLogical,
                      G4LogicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
                const G4double width,
                const G4double offset = 0.);

    /**
     * Similar to the constructor above, except for the mother pointer's type
     * being here a G4VPhysicalVolume.
     *  @param[in] pName The volume name.
     *  @param[in] pLogical Pointer to the logical volume of the replica.
     *  @param[in] pMother Pointer to the physical volume of the mother.
     *  @param[in] pAxis The axis along which do the replication.
     *  @param[in] nReplicas The number of copies to replicate.
     *  @param[in] width The witdh of the replicated object along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     */
    G4PVReplica(const G4String& pName,
                      G4LogicalVolume* pLogical,
                      G4VPhysicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
                const G4double width,
                const G4double offset = 0.);

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4PVReplica(__void__&);

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4PVReplica(const G4PVReplica&) =  delete;
    G4PVReplica& operator=(const G4PVReplica&) = delete;

    /**
     * Default Destructor.
     */
    ~G4PVReplica() override = default;

    /**
     * Returns the volume type characterisation.
     */
    EVolume VolumeType() const override;

    /**
     * Not used.
     */
    G4bool IsMany() const override;

    /**
     * Returns true.
     */
    G4bool IsReplicated() const override;

    /**
     * Accessor/modifier for optional handling of the copy-number.
     */
    G4int GetCopyNo() const override;
    void  SetCopyNo(G4int CopyNo) override;

    /**
     * Returns false and nullptr.
     */
    G4bool IsParameterised() const override;
    G4VPVParameterisation* GetParameterisation() const override;

    /**
     * Returns the number of replications.
     */
    G4int GetMultiplicity() const override;

    /**
     * Fills arguments with the attributes from the base replica.
     * @param[in,out] axis Axis of parameterisation returned.
     * @param[in,out] nReplicas The number of replica copies.
     * @param[in,out] width Width of the replica object.
     * @param[in,out] offset Potential offset in replication.
     * @param[in,out] consuming Flag of replica characterisation (always true
     *                for pure replicas).
     */
    void GetReplicationData(EAxis& axis,
                            G4int& nReplicas,
                            G4double& width,
                            G4double& offset,
                            G4bool& consuming) const override;

    /**
     * Sets a unique code for each type of regular structure.
     *  @note It must be called only during detector construction.
     *        It can also be used to prepare any corresponding special
     *        navigation 'conditions'.
     */
    virtual void SetRegularStructureId( G4int code ); 

    /**
     * Accessors for specialised geometries.
     */
    G4bool IsRegularStructure() const override; 
    G4int GetRegularStructureId() const override;

    // Methods for handling of MT instances

    /**
     * Returns the MT instance ID.
     */
    inline G4int GetInstanceID() const  { return instanceID; }

    /**
     * Returns the private data instance manager.
     */
    static const G4PVRManager& GetSubInstanceManager();

    /**
     * This method is similar to the constructor. It is used by each worker
     * thread to achieve the partial effect as that of the master thread.
     */
    void InitialiseWorker(G4PVReplica* pMasterObject);

    /**
     * This method is similar to the destructor. It is used by each worker
     * thread to achieve the partial effect as that of the master thread.
     */
    void TerminateWorker(G4PVReplica* pMasterObject);

  protected:

    /**
     * Constructor for derived type(s): PVParameterised, PVDivision, ...
     * Does not set mother or register in mother volume -- leaves it to
     * derived type.
     */
    G4PVReplica(const G4String& pName,
                      G4int nReplicas,
                      EAxis pAxis,
                      G4LogicalVolume*  pLogical,
                      G4LogicalVolume*  pMotherLogical);

  private:

    /**
     * Performs sanity checks on parameters and initialises data.
     */
    void CheckAndSetParameters(const EAxis pAxis, const G4int nReplicas,
                               const G4double width, const G4double offset);

    /**
     * Checks that this volume is the only daughter of its proposed mother
     * volume.
     */
    void CheckOnlyDaughter(G4LogicalVolume* pMotherLogical);
   
  protected:

    EAxis faxis;
    G4int fnReplicas;
    G4double fwidth, foffset;
 
  private:

    G4int fRegularVolsId = 0;

    G4int instanceID;  // Used as instance ID
    G4GEOM_DLL static G4PVRManager subInstanceManager; // Uses G4PVRManager
};

/**
 * @note:
 *
 * The type G4PVRManager is introduced to encapsulate the methods used by
 * both the master thread and worker threads to allocate memory space for
 * the fields encapsulated by the class G4ReplicaData. When each thread
 * initializes the value for these fields, it refers to them using a macro
 * definition defined below. For every G4PVReplica instance, there is
 * a corresponding G4ReplicaData instance. All G4ReplicaData instances are
 * organized by the class G4PVRManager as an array.
 * The field "int instanceID" is added to the class G4PVReplica.
 * The value of this field in each G4LogicalVolume instance is the subscript
 * of the corresponding G4ReplicaData instance.
 * In order to use the class  G4PVRManager, we add a static member in the
 * class G4LogicalVolume as follows: "static G4PVRManager subInstanceManager".
 * For the master thread, the array for G4ReplicaData instances grows
 * dynamically along with G4PVReplica instances arecreated.
 * For each worker thread, it copies the array of G4ReplicaData instances
 * from the master thread.
 * In addition, it invokes a method similiar to the constructor explicitly
 * to achieve the partial effect for each instance in the array.
 */

#endif
