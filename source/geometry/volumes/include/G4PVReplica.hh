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
// $Id: G4PVReplica.hh 85846 2014-11-05 15:45:28Z gcosmo $
//
// 
// class G4PVReplica
//
// Class description:
//
// Represents many touchable detector elements differing only in their
// positioning. The elements' positions are calculated by means of a simple
// linear formula, and the elements completely fill the containing mother
// volume.
// 
// G4PVReplica(const G4String& pName,
//                   G4LogicalVolume *pLogical,
//                   G4LogicalVolume *pMother,
//             const EAxis pAxis,
//             const G4int nReplicas,
//             const G4double width,
//             const G4double offset=0);
//
// Replication may occur along:
//
// o Cartesian axes (kXAxis,kYAxis,kZAxis)
//
//   The replications, of specified width have coordinates of
//   form (-width*(nReplicas-1)*0.5+n*width,0,0) where n=0.. nReplicas-1
//   for the case of kXAxis, and are unrotated.
//
// o Radial axis (cylindrical polar) (kRho)
//
//   The replications are cons/tubs sections, centred on the origin
//   and are unrotated.
//   They have radii of width*n+offset to width*(n+1)+offset
//                      where n=0..nReplicas-1
//
// o Phi axis (cylindrical polar) (kPhi)
//   The replications are `phi sections' or wedges, and of cons/tubs form
//   They have phi of offset+n*width to offset+(n+1)*width where
//   n=0..nReplicas-1

// History:
// 29.07.95 P.Kent           - First non-stub version
// 26.10.97 J.Apostolakis    - Added constructor that takes mother LV
// 16.02.98 J.Apostolakis    - Added copy number
// 13.01.13 G.Cosmo, A.Dotti - Modified for thread-safety for MT
// ----------------------------------------------------------------------
#ifndef G4PVREPLICA_HH
#define G4PVREPLICA_HH

#include "G4VPhysicalVolume.hh"
#include "G4GeomSplitter.hh"

class G4ReplicaData
{
  // Encapsulates the fields of the class G4PVReplica that may not be
  // read-only. G4PVReplica inherits from the class G4VPhysicalVolume.
  // The fields from the ancestor that may not be read-only are handled
  // by the ancestor class.

public:

  void initialize() {}

  G4int    fcopyNo;
};

// The type G4PVRManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4ReplicaData. When each thread
// initializes the value for these fields, it refers to them using a macro
// definition defined below. For every G4PVReplica instance, there is
// a corresponding G4ReplicaData instance. All G4ReplicaData instances are
// organized by the class G4PVRManager as an array.
// The field "int instanceID" is added to the class G4PVReplica.
// The value of this field in each G4LogicalVolume instance is the subscript
// of the corresponding G4ReplicaData instance.
// In order to use the class  G4PVRManager, we add a static member in the
// class G4LogicalVolume as follows: "static G4PVRManager subInstanceManager".
// For the master thread, the array for G4ReplicaData instances grows
// dynamically along with G4PVReplica instances arecreated.
// For each worker thread, it copies the array of G4ReplicaData instances
// from the master thread.
// In addition, it invokes a method similiar to the constructor explicitly
// to achieve the partial effect for each instance in the array.
//
typedef G4GeomSplitter<G4ReplicaData> G4PVRManager;

// This macro changes the references to fields that are now encapsulated
// in the class G4ReplicaData.
//
#define G4MT_copyNo ((subInstanceManager.offset[instanceID]).fcopyNo)

class G4PVReplica : public G4VPhysicalVolume
{
  public:  // with description

    G4PVReplica(const G4String& pName,
                      G4LogicalVolume* pLogical,
                      G4LogicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
                const G4double width,
                const G4double offset=0);

  public:  // without description

    G4PVReplica(const G4String& pName,
                      G4LogicalVolume* pLogical,
                      G4VPhysicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
                const G4double width,
                const G4double offset=0);

    G4PVReplica(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  public:  // with description

    virtual ~G4PVReplica();

    G4bool IsMany() const;
    G4bool IsReplicated() const;

    virtual G4int GetCopyNo() const;
    virtual void  SetCopyNo(G4int CopyNo);
    virtual G4bool IsParameterised() const;
    virtual G4VPVParameterisation* GetParameterisation() const;
    virtual G4int GetMultiplicity() const;
    virtual void GetReplicationData(EAxis& axis,
                                    G4int& nReplicas,
                                    G4double& width,
                                    G4double& offset,
                                    G4bool& consuming) const;

    virtual void SetRegularStructureId( G4int Code ); 
      // This method must set a unique code for each type of regular structure.
      // - It must be called only during detector construction.
      // - It can also be used to prepare any corresponding special
      //  navigation 'conditions'.

    G4bool IsRegularStructure() const; 
    G4int GetRegularStructureId() const;
      // Accessors for specialised geometries

  public:  // without description

    inline G4int GetInstanceID() const  { return instanceID; }
      // Returns the instance ID.

    static const G4PVRManager& GetSubInstanceManager();
      // Returns the private data instance manager.

    void InitialiseWorker(G4PVReplica *pMasterObject);
      // This method is similar to the constructor. It is used by each worker
      // thread to achieve the partial effect as that of the master thread.

    void TerminateWorker(G4PVReplica *pMasterObject);
      // This method is similar to the destructor. It is used by each worker
      // thread to achieve the partial effect as that of the master thread.

  private:

    void CheckAndSetParameters(const EAxis pAxis, const G4int nReplicas,
                               const G4double width, const G4double offset);
    G4PVReplica(const G4PVReplica&);
    G4PVReplica& operator=(const G4PVReplica&);

  protected:

    EAxis faxis;
    G4int fnReplicas;
    G4double fwidth,foffset;
 
  private:

    G4int fRegularStructureCode; 
    G4int fRegularVolsId;

    G4int instanceID;
      // This new field is used as instance ID.
    G4GEOM_DLL static G4PVRManager subInstanceManager;
      // This new field helps to use the class G4PVRManager introduced above.
};

#endif
