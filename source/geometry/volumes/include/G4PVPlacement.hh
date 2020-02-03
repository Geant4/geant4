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
// G4PVPlacement
//
// Class description:
//
// Class representing a single volume positioned within and relative
// to a mother volume.

// 24.07.95 P.Kent, First non-stub version
// 25.07.96 P.Kent, Modified interface for new `Replica' capable geometry
// 28.08.96 P.Kent. Tidied + transform replaced by rotmat+vector
// ----------------------------------------------------------------------
#ifndef G4PVPLACEMENT_HH
#define G4PVPLACEMENT_HH

#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"

class G4PVPlacement : public G4VPhysicalVolume
{
  public:  // with description

    G4PVPlacement(G4RotationMatrix* pRot,
            const G4ThreeVector& tlate,
                  G4LogicalVolume* pCurrentLogical,
            const G4String& pName,
                  G4LogicalVolume* pMotherLogical,
                  G4bool pMany,
                  G4int  pCopyNo,
                  G4bool pSurfChk = false);
      // Initialise a single volume, positioned in a frame which is rotated by
      // *pRot and traslated by tlate, relative to the coordinate system of the
      // mother volume pMotherLogical.
      // If pRot=nullptr the volume is unrotated with respect to its mother.
      // The physical volume is added to the mother's logical volume.
      // Arguments particular to G4PVPlacement:
      //   pMany Currently NOT used. For future use to identify if the volume
      //         is meant to be considered an overlapping structure, or not.
      //   pCopyNo should be set to 0 for the first volume of a given type.
      //   pSurfChk if true activates check for overlaps with existing volumes.
      // This is a very natural way of defining a physical volume, and is
      // especially useful when creating subdetectors: the mother volumes are
      // not placed until a later stage of the assembly program.

    G4PVPlacement(const G4Transform3D& Transform3D,
                        G4LogicalVolume* pCurrentLogical,
                  const G4String& pName,
                        G4LogicalVolume* pMotherLogical,
                        G4bool pMany,
                        G4int pCopyNo,
                        G4bool pSurfChk = false);
      // Additional constructor, which expects a G4Transform3D that represents 
      // the direct rotation and translation of the solid (NOT of the frame).  
      // The G4Transform3D argument should be constructed by:
      //  i) First rotating it to align the solid to the system of 
      //     reference of its mother volume *pMotherLogical, and 
      // ii) Then placing the solid at the location Transform3D.getTranslation(),
      //     with respect to the origin of the system of coordinates of the
      //     mother volume.  
      // [ This is useful for the people who prefer to think in terms 
      //   of moving objects in a given reference frame. ]
      // All other arguments are the same as for the previous constructor.

  public:  // without description

    G4PVPlacement(G4RotationMatrix* pRot,
            const G4ThreeVector& tlate,
            const G4String& pName,
                  G4LogicalVolume* pLogical,
                  G4VPhysicalVolume* pMother,
                  G4bool pMany,
                  G4int pCopyNo,
                  G4bool pSurfChk = false);
      // A simple variation of the 1st constructor, only specifying the
      // mother volume as a pointer to its physical volume instead of its
      // logical volume. The effect is exactly the same.

    G4PVPlacement(const G4Transform3D& Transform3D,
                  const G4String& pName,
                        G4LogicalVolume* pLogical,
                        G4VPhysicalVolume* pMother,
                        G4bool pMany,
                        G4int pCopyNo,
                        G4bool pSurfChk = false);
      // Utilises both variations above (from 2nd and 3rd constructor).
      // The effect is the same as for the 2nd constructor.

  public:  // with description

    virtual ~G4PVPlacement();
      // Default destructor.

    inline G4int GetCopyNo() const  { return fcopyNo; }

    void SetCopyNo(G4int CopyNo);
      // Gets and sets the copy number of the volume.

    G4bool CheckOverlaps(G4int res = 1000, G4double tol = 0.,
                         G4bool verbose = true, G4int maxErr = 1);
      // Verifies if the placed volume is overlapping with existing
      // daughters or with the mother volume. Provides default resolution
      // for the number of points to be generated and verified.
      // A tolerance for the precision of the overlap check can be specified,
      // by default it is set to maximum precision.
      // Reports a maximum of overlaps errors according to parameter in input.
      // Returns true if the volume is overlapping.

  public:  // without description

    G4PVPlacement(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4PVPlacement(const G4PVPlacement&) = delete;
    G4PVPlacement& operator=(const G4PVPlacement&) = delete;
      // Copy constructor and assignment operator not allowed.

    G4bool IsMany() const;
    G4bool IsReplicated() const;
    G4bool IsParameterised() const;
    G4VPVParameterisation* GetParameterisation() const;
    void GetReplicationData(EAxis& axis,
                            G4int& nReplicas,
                            G4double& width,
                            G4double& offset,
                            G4bool& consuming) const;
    G4bool IsRegularStructure() const; 
    G4int GetRegularStructureId() const; 

    EVolume VolumeType() const;
   
  private:

    static G4RotationMatrix* NewPtrRotMatrix(const G4RotationMatrix& RotMat);
      // Auxiliary function for 2nd constructor (one with G4Transform3D).
      // Creates a new RotMatrix on the heap (using "new") and copies 
      // its argument into it.

  private:

    G4bool fmany = false;   // flag for overlapping structure - not used
    G4bool fallocatedRotM = false;  // flag for allocation of Rotation Matrix
    G4int fcopyNo = 0;      // for identification

};

#endif

