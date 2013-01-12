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
// $Id$
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
// 29.07.95 P.Kent         First non-stub version
// 26.10.97 J.Apostolakis  Added constructor that takes mother logical volume
// 16.02.98 J.Apostolakis  Added copy number
// ----------------------------------------------------------------------
#ifndef G4PVREPLICA_HH
#define G4PVREPLICA_HH

#include "G4VPhysicalVolume.hh"

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The class ReplicaPrivateSubclass is introduced to encapsulate
//the fields of the class G4PVReplica that may not be read-only.
//G4PVReplica inherits from the class G4VPhysicalVolume. The fields
//from the ancestor that may not be read-only are handled the ancestor
//class.
#ifndef REPLICASHAREDSUBCLASS_HH
#define REPLICASHAREDSUBCLASS_HH

class ReplicaPrivateSubclass
{
public:
  G4int    fcopyNo;
  void initialize(){};
};
#endif

//01.25.2009 Xin Dong: Phase II change for Geant4 multithreading.
//The class G4PVReplicaSubInstanceManager is introduced to
//encapsulate the methods used by both the master thread and
//worker threads to allocate memory space for the fields encapsulated
//by the class ReplicaPrivateSubclass. When each thread
//initializes the value for these fields, it refers to them using a macro
//definition defined below. For every G4PVReplica instance, there is
//a corresponding ReplicaPrivateSubclass instance. All
//ReplicaPrivateSubclass instances are organized by the
//class G4PVReplicaSubInstanceManager as an array. The field "
//int g4pvreplicaInstanceID" is added to the class G4PVReplica.
//The value of this field in each G4LogicalVolume instance is the subscript
//of the corresponding ReplicaPrivateSubclass instance. In order
//to use the class  G4PVReplicaSubInstanceManager, we add a static member in
//the class G4LogicalVolume as follows: "
//static G4PVReplicaSubInstanceManager g4pvreplicaSubInstanceManager".
//For the master thread, the array for ReplicaPrivateSubclass
//instances grows dynamically along with G4PVReplica instances are
//created. For each worker thread, it copies the array of
//ReplicaPrivateSubclass instances from the master thread.
//In addition, it invokes a method similiar to the constructor explicitly
//to achieve the partial effect for each instance in the array.
#ifndef G4PVREPLICASUBINSTANCEMANAGER_HH
#define G4PVREPLICASUBINSTANCEMANAGER_HH

#include "G4MTTransitory.hh"

typedef G4MTPrivateSubInstanceManager<ReplicaPrivateSubclass> G4PVReplicaSubInstanceManager;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//These macros changes the references to fields that are now encapsulated
//in the class ReplicaPrivateSubclass.
#define fcopyNoG4MTThreadPrivate ((g4pvreplicaSubInstanceManager.offset[g4pvreplicaInstanceID]).fcopyNo)

#endif


class G4PVReplica : public G4VPhysicalVolume
{

  public:

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This new field is used as instance ID.
    int g4pvreplicaInstanceID;

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This new field helps to use the class G4PVReplicaSubInstanceManager
    //introduced above.
    static G4PVReplicaSubInstanceManager g4pvreplicaSubInstanceManager;

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This method is similar to the constructor. It is used by each worker
    //thread to achieve the partial effect as that of the master thread.
    void SlaveG4PVReplica(G4PVReplica *pMasterObject);

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This method is similar to the destructor. It is used by each worker
    //thread to achieve the partial effect as that of the master thread.
    void DestroySlaveG4PVReplica(G4PVReplica *pMasterObject);

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

  protected:

    EAxis faxis;
    G4int fnReplicas;
    G4double fwidth,foffset;
 
    // 01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    // This field is move from the original class definition to be
    // encapsulated by the class ReplicaPrivateSubclass.
    //    G4int    fcopyNoG4MTThreadPrivate;

  private:

    void CheckAndSetParameters(const EAxis pAxis, const G4int nReplicas,
                               const G4double width, const G4double offset);
    G4PVReplica(const G4PVReplica&);
    const G4PVReplica& operator=(const G4PVReplica&);

  private:

    G4int fRegularStructureCode; 
    G4int fRegularVolsId;
};

#endif
