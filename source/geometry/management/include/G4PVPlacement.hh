// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVPlacement.hh,v 1.3 2000-04-20 16:49:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVPlacement
//
// Class description:
//
// Class representing a single volume positioned within and relative
// to a mother volume.

// History:
// 24.07.95 P.Kent   First non-stub version
// 25.07.96 P.Kent   Modified interface for new `Replica' capable geometry
// 28.08.96 P.Kent   Tidied + transform replaced by rotmat+vector
// 28.02.97 J.Apostolakis Added 2nd constructor with G4Transform3D of solid.
// 11.07.97 J.Apostolakis Added 3rd constructor with pMotherLogical 
// 11.05.98 J.Apostolakis Added 4th constructor with G4Transform3D & pMotherLV

#ifndef G4PVPLACEMENT_HH
#define G4PVPLACEMENT_HH

#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"

class G4PVPlacement : public G4VPhysicalVolume
{
  public:  // with description

    G4PVPlacement(G4RotationMatrix *pRot,
		  const G4ThreeVector &tlate,
		  const G4String &pName,
		  G4LogicalVolume *pLogical,
		  G4VPhysicalVolume *pMother,
		  G4bool pMany,
		  G4int pCopyNo);
      // Initialise a single volume, positioned in a frame which is rotated by
      // *pFrameRot, relative to the coordinate system of the mother volume
      // pMother. The center of the object is then placed at volumeCenterCrd
      // in the new coordinates. 
      // If pRot=0 the volume is unrotated with respect to its mother.
      // The physical volume is added to the mother's logical volume.
      //   (The above are exactly the arguments of G4VPhysicalVolume)
      // Arguments particular to G4PVPlacement:
      //   pMany Currently NOT used. For future use to identify if the volume
      //         is MANY in the GEANT 3 sense, or not.
      //   pCopyNo should be set to 0 for the first volume of a given type.

    G4PVPlacement(const G4Transform3D &Transform3D,
		  const G4String &pName,
		  G4LogicalVolume *pLogical,
		  G4VPhysicalVolume *pMother,
		  G4bool pMany,
		  G4int pCopyNo);
      // Additional constructor, which expects a G4Transform3D that represents 
      // the direct rotation and translation of the solid (NOT of the frame).  
      // The G4Transform3D argument should be constructed by:
      //  i) First rotating it to align the solid to the system of 
      //      reference of its mother volume *pMother, and 
      // ii) Then placing the solid at the location Transform3D.getTranslation(),
      //      with respect to the origin of the system of coordinates of the
      //      mother volume.  
      // [ This is useful for the people who prefer to think in terms 
      // of moving objects in a given reference frame. ]
      // All other arguments are the same as for the previous constructor.

    G4PVPlacement(G4RotationMatrix *pRot,
		  const G4ThreeVector &tlate,
		  G4LogicalVolume *pCurrentLogical,
		  const G4String& pName,
		  G4LogicalVolume *pMotherLogical,
		  G4bool pMany,
		  G4int pCopyNo);
      // A simple variation of the 1st constructor, only specifying the
      // mother volume as a pointer to its logical volume instead of its
      // physical volume. [ This is a very natural way of defining a physical
      // volume, and is especially useful when creating subdetectors: the
      // mother volumes is not placed until a later stage of the assembly
      // program. ]

    G4PVPlacement(const G4Transform3D &Transform3D,
                  G4LogicalVolume *pCurrentLogical,
                  const G4String& pName,
                  G4LogicalVolume *pMotherLogical,
                  G4bool pMany,
                  G4int pCopyNo);
      // Utilises both variations above (from 2nd and 3rd constructor).

    ~G4PVPlacement();
      // Default destructor.

    virtual G4int GetCopyNo() const;
    virtual void  SetCopyNo(G4int CopyNo);
      // Gets and sets the copy number of the volume.

    virtual void Setup(G4VPhysicalVolume *pMother);
      // Sets the pointer to the mother volume.
      // Must not be called when geometry closed.

  public:  // without description

    virtual G4bool IsMany() const;
    virtual G4bool IsReplicated() const;
    virtual G4VPVParameterisation* GetParameterisation() const;
    virtual void GetReplicationData(EAxis& axis,
                                   G4int& nReplicas,
				   G4double& width,
                                   G4double& offset,
                                   G4bool& consuming) const;

  private:

    G4bool fmany;	    // flag for booleans/MANY - not used
    G4bool fallocatedRotM;  // flag for allocation of Rotation Matrix
    G4int fcopyNo;	    // for identification

    static G4RotationMatrix* NewPtrRotMatrix(const G4RotationMatrix &RotMat);
      // Auxiliary function for 2nd constructor (one with G4Transform3D).
      // Creates a new RotMatrix on the heap (using "new") and copies 
      // its argument into it.
};

#endif

