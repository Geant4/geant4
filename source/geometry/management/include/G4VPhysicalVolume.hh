// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPhysicalVolume.hh,v 1.4 2000-04-20 16:49:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifndef G4VPHYSICALVOLUME_HH
#define G4VPHYSICALVOLUME_HH

#include "globals.hh"
#include "geomdefs.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPVParameterisation;

class G4VPhysicalVolume
{
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

    G4bool operator == (const G4VPhysicalVolume& p) const;
      // Define equality by equal addresses only.

    // Access functions

      // The following are accessor functions that make a distinction
      // between whether the rotation/translation is being made for the
      // frame or the object/volume that is being placed.
      // (They are the inverse of each other).
    G4RotationMatrix* GetObjectRotation() const;       //  Obsolete 
    G4RotationMatrix  GetObjectRotationValue() const;  //  Replacement
    G4ThreeVector  GetObjectTranslation() const;
      // Return the rotation/translation of the Object relative to the mother.
    const G4RotationMatrix* GetFrameRotation() const;
    G4ThreeVector  GetFrameTranslation() const;
      // Return the rotation/translation of the Frame used to position 
      // this volume in its mother volume (opposite of object rot/trans).

    // Older access functions, that do not distinguish between frame/object!

    const G4ThreeVector& GetTranslation() const;
    const G4RotationMatrix* GetRotation() const;
      // Old access functions, that do not distinguish between frame/object!
      // They return the translation/rotation of the volume.

    // Set functions

    void SetTranslation(const G4ThreeVector &v);
    G4RotationMatrix* GetRotation();
    void SetRotation(G4RotationMatrix*);
      // NOT INTENDED FOR GENERAL USE.
      // Non constant versions of above. Used to change transformation
      // for replication/parameterisation mechanism.

    G4LogicalVolume* GetLogicalVolume() const;
      // Return the associated logical volume.
    void SetLogicalVolume(G4LogicalVolume *pLogical);
      // Set the logical volume. Must not be called when geometry closed.

    G4VPhysicalVolume* GetMother() const;
      // Return the current mother pointer.
    void SetMother(G4VPhysicalVolume *pMother);
      // Set the mother volume. Must not be called when geometry closed.

    G4String GetName() const;
      // Return the volume's name.
    void SetName(const G4String& pName);
      // Set the volume's name.

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
    virtual G4VPVParameterisation* GetParameterisation() const = 0;
      // Return replicas parameterisation object (able to compute dimensions
      // and transformations of replicas), or NULL if not applicable.
    virtual void GetReplicationData(EAxis& axis,
                                   G4int& nReplicas,
				   G4double& width,
                                   G4double& offset,
                                   G4bool& consuming) const = 0;
      // Return replication information. No-op for no replicated volumes.
    virtual void Setup(G4VPhysicalVolume *pMother) = 0;
      // Perform any initialisation/setup necessary for the given volume.
      // [Set the current mother pointer to refer to the specified mother,
      // by calling SetMother]

  protected:

    G4RotationMatrix *frot;
    G4ThreeVector ftrans;

  private:

    G4LogicalVolume *flogical;	// The logical volume representing the
				// physical and tracking attributes of
                                // the volume
    G4String fname;		// name of the volume
    G4VPhysicalVolume *fmother;	// The current moher volume
};

#include "G4VPhysicalVolume.icc"

#endif
