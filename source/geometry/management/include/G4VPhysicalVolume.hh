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
// class G4VPhysicalVolume
//
// Class description:
//
// This is an Abstract Base class for the representation of positioned volume.  
// The volume is placed within a mother volume,  relative to its coordinate 
// system.  Either a single positioned volume or many positioned volume can 
// be represented by a particular G4VPhysicalVolume.

// History:
// 09.11.99 J.Apostolakis  Added GetObjectRotationValue() method & redid comments.
// 28.08.96 P.Kent Replaced transform by rotmat + vector
// 25.07.96 P.Kent Modified interface for new `Replica' capable geometry 
// 24.07.95 P.Kent First non-stub version
// --------------------------------------------------------------------
#ifndef G4VPHYSICALVOLUME_HH
#define G4VPHYSICALVOLUME_HH

#include "G4Types.hh"
#include "G4String.hh"

#include "geomdefs.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPVParameterisation;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The class PhysicalVolumePrivateSubclass is introduced to
//encapsulate the fields of the class G4VPhysicalVolume that may not
//be read-only.          
#ifndef PHYSICALVOLUMEPRIVATESUBCLASS_HH
#define PHYSICALVOLUMEPRIVATESUBCLASS_HH

class PhysicalVolumePrivateSubclass
{
public:
  G4RotationMatrix *frot;
  G4ThreeVector ftrans;
  void initialize() {};
};
#endif

//01.25.2009 Xin Dong: Phase II change for Geant4 multithreading.
//The class G4VPhysicalVolumeSubInstanceManager is introduced to
//encapsulate the methods used by both the master thread and
//worker threads to allocate memory space for the fields encapsulated
//by the class PhysicalVolumePrivateSubclass. When each thread
//initializes the value for these fields, it refers to them using a macro
//definition defined below. For every G4VPhysicalVolume instance, there is
//a corresponding PhysicalVolumePrivateSubclass instance. All
//PhysicalVolumePrivateSubclass instances are organized by the
//class G4VPhysicalVolumeSubInstanceManager as an array. The field "
//int g4vphysicalVolumeInstanceID" is added to the class G4VPhysicalVolume.
//The value of this field in each G4VPhysicalVolume instance is the subscript
//of the corresponding PhysicalVolumePrivateSubclass instance. In order
//to use the class G4VPhysicalVolumeSubInstanceManager, we add a static member in      
//the class G4VPhysicalVolume as follows: "                                       
//static G4VPhysicalVolumeSubInstanceManager g4vphysicalVolumeSubInstanceManager;".
//For the master thread, the array for PhysicalVolumePrivateSubclass
//instances grows dynamically along with G4VPhysicalVolume instances are
//created. For each worker thread, it copies the array of
//PhysicalVolumePrivateSubclass instances from the master thread.           
//In addition, it invokes a method similiar to the constructor explicitly       
//to achieve the partial effect for each instance in the array. 
#ifndef G4VPHYSICALVOLUMESUBINSTANCEMANAGER_HH
#define G4VPHYSICALVOLUMESUBINSTANCEMANAGER_HH

#include "G4MTTransitory.hh"
typedef G4MTPrivateSubInstanceManager<PhysicalVolumePrivateSubclass> G4VPhysicalVolumeSubInstanceManager;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//These macros changes the references to fields that are now encapsulated
//in the class PhysicalVolumePrivateSubclass.
#define frotG4MTThreadPrivate ((g4vphysicalVolumeSubInstanceManager.offset[g4vphysicalVolumeInstanceID]).frot)
#define ftransG4MTThreadPrivate ((g4vphysicalVolumeSubInstanceManager.offset[g4vphysicalVolumeInstanceID]).ftrans)

#endif

class G4VPhysicalVolume
{
  public:

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This new field is used as instance ID.
    int g4vphysicalVolumeInstanceID;

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This new field helps to use the class G4VPhysicalVolumeSubInstanceManager
    //introduced above.
    static G4VPhysicalVolumeSubInstanceManager g4vphysicalVolumeSubInstanceManager;

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This method is similar to the constructor. It is used by each worker
    //thread to achieve the partial effect as that of the master thread.
    void SlaveG4VPhysicalVolume(G4VPhysicalVolume *pMasterObject, G4RotationMatrix *pRot,
				const G4ThreeVector &tlate);

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
    //This method is similar to the destructor. It is used by each worker
    //thread to achieve the partial effect as that of the master thread.
    void DestroySlaveG4VPhysicalVolume(G4VPhysicalVolume *pMasterObject);

  public:  // with description

    G4VPhysicalVolume(G4RotationMatrix *pRot,
                const G4ThreeVector &tlate,
                const G4String &pName,
                      G4LogicalVolume *pLogical,
                      G4VPhysicalVolume *pMother);


      // Initialise volume, positioned in a frame which is rotated by *pRot, 
      // relative to the coordinate system of the mother volume pMother.
      // The center of the object is then placed at tlate in the new
      // coordinates. If pRot=0 the volume is unrotated with respect to its
      // mother. The physical volume is added to the mother's logical volume.
      //
      // Must be called by all subclasses. pMother must point to a valid parent
      // volume, except in the case of the world/top volume, when it =0.
      // 
      // The constructor also registers volume with physical volume Store.
      // Note that the Store may be removed or dynamically built in future
      // because of memory constraints.

    virtual ~G4VPhysicalVolume();
      // Destructor, will be subclassed. Removes volume from volume Store.

    inline G4bool operator == (const G4VPhysicalVolume& p) const;
      // Equality defined by equal addresses only.

    // Access functions

      // The following are accessor functions that make a distinction
      // between whether the rotation/translation is being made for the
      // frame or the object/volume that is being placed.
      // (They are the inverse of each other).
    G4RotationMatrix* GetObjectRotation() const;              //  Obsolete 
    inline G4RotationMatrix  GetObjectRotationValue() const;  //  Replacement
    inline G4ThreeVector  GetObjectTranslation() const;
      // Return the rotation/translation of the Object relative to the mother.
    inline const G4RotationMatrix* GetFrameRotation() const;
    inline G4ThreeVector  GetFrameTranslation() const;
      // Return the rotation/translation of the Frame used to position 
      // this volume in its mother volume (opposite of object rot/trans).

    // Older access functions, that do not distinguish between frame/object!

    inline const G4ThreeVector& GetTranslation() const;
    inline const G4RotationMatrix* GetRotation() const;
      // Old access functions, that do not distinguish between frame/object!
      // They return the translation/rotation of the volume.

    // Set functions

    inline void SetTranslation(const G4ThreeVector &v);
    inline G4RotationMatrix* GetRotation();
    inline void SetRotation(G4RotationMatrix*);
      // NOT INTENDED FOR GENERAL USE.
      // Non constant versions of above. Used to change transformation
      // for replication/parameterisation mechanism.

    inline G4LogicalVolume* GetLogicalVolume() const;
      // Return the associated logical volume.
    inline void SetLogicalVolume(G4LogicalVolume *pLogical);
      // Set the logical volume. Must not be called when geometry closed.

    inline G4LogicalVolume* GetMotherLogical() const;
      // Return the current mother logical volume pointer.
    inline void SetMotherLogical(G4LogicalVolume *pMother);
      // Set the mother logical volume. Must not be called when geometry closed.

    inline const G4String& GetName() const;
      // Return the volume's name.
    inline void SetName(const G4String& pName);
      // Set the volume's name.

    virtual G4int GetMultiplicity() const;
      // Returns number of object entities (1 for normal placements,
      // n for replicas or parameterised).

    // Functions required of subclasses

    virtual G4bool IsMany() const = 0;
      // Return true if the volume is MANY (not implemented yet).
    virtual G4int GetCopyNo() const = 0;
      // Return the volumes copy number.
    virtual void  SetCopyNo(G4int CopyNo) = 0;
      // Set the volumes copy number.
    virtual G4bool IsReplicated() const = 0;
      // Return true if replicated (single object instance represents
      // many real volumes), else false.
    virtual G4bool IsParameterised() const = 0;
      // Return true if parameterised (single object instance represents
      // many real parameterised volumes), else false.
    virtual G4VPVParameterisation* GetParameterisation() const = 0;
      // Return replicas parameterisation object (able to compute dimensions
      // and transformations of replicas), or NULL if not applicable.
    virtual void GetReplicationData(EAxis& axis,
                                    G4int& nReplicas,
                                    G4double& width,
                                    G4double& offset,
                                    G4bool& consuming) const = 0;
      // Return replication information. No-op for no replicated volumes.
    virtual G4bool  IsRegularStructure() const = 0;
      // Returns true if the underlying volume structure is regular.
    virtual G4int  GetRegularStructureId() const = 0;
      // Returns non-zero code in case the underlying volume structure 
      //  is regular, voxel-like.  Value is id for structure type.
      //  If non-zero the volume is a candidate for specialised 
      //  navigation such as 'nearest neighbour' directly on volumes.
    virtual G4bool CheckOverlaps(G4int res=1000, G4double tol=0.,
                                 G4bool verbose=true);
      // Verifies if the placed volume is overlapping with existing
      // daughters or with the mother volume. Provides default resolution
      // for the number of points to be generated and verified.
      // Concrete implementation is done and required only for placed and
      // parameterised volumes. Returns true if the volume is overlapping.

  public:  // without description

    G4VPhysicalVolume(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  private:

    G4VPhysicalVolume(const G4VPhysicalVolume&);
    G4VPhysicalVolume& operator=(const G4VPhysicalVolume&);
      // Private copy constructor and assignment operator.

  protected:

  // 01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  // This field is move from the original class definition to be
  // encapsulated by the class PhysicalVolumePrivateSubclass.
  //    G4RotationMatrix *frotG4MTThreadPrivate;

  // 01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
  // This field is move from the original class definition to be
  // encapsulated by the class PhysicalVolumePrivateSubclass.
  //    G4ThreeVector ftransG4MTThreadPrivate;

  private:

    G4LogicalVolume *flogical;   // The logical volume representing the
                                 // physical and tracking attributes of
                                 // the volume
    G4String fname;              // The name of the volume
    G4LogicalVolume   *flmother; // The current mother logical volume
};

#include "G4VPhysicalVolume.icc"

#endif
