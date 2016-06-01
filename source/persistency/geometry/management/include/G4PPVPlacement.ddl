// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPVPlacement.ddl,v 2.3 1998/11/10 18:29:53 morita Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
//          P-versieon of G4PVPlacement      Takashi.Sasaki@kek.jp
#ifndef G4PPVPLACEMENT_DDL
#define G4PPVPLACEMENT_DDL 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PLogicalVolume;

#include "G4PVPhysicalVolume.hh"

class G4PPVPlacement : public G4PVPhysicalVolume
{
public:
  G4PPVPlacement( G4VPhysicalVolume *PhysVol,
		          HepRef(G4PLogicalVolume) persLogVol);

  ~G4PPVPlacement();

  G4VPhysicalVolume* MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother );

  virtual G4bool IsMany() const;
  virtual G4int GetCopyNo() const;

protected:
  virtual void  SetCopyNo(G4int CopyNo);

private:
  G4bool fmany;	    // flag for booleans
//  G4bool fallocatedRotM;  // flag for allocation of Rotation Matrix
  G4int fcopyNo;	    // for identification
};

#endif
