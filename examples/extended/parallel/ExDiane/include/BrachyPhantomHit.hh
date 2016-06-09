//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: BrachyPhantomHit.hh,v 1.2 2004/05/25 08:36:17 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
//    ********************************
//    *                              *
//    *    BrachyPhantomHit.hh       *
//    *                              *
//    ********************************
//
//Code developed by: S. Guatelli
//
// It manages the hits and the enrgy deposit in the phantom associated with 
// each hit ...

#ifndef BrachyPhantomHit_h
#define BrachyPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class BrachyPhantomHit : public G4VHit
{
public:
  BrachyPhantomHit(G4LogicalVolume* ,G4int ,G4int ,G4int );
  ~BrachyPhantomHit();
  BrachyPhantomHit(const BrachyPhantomHit &right);
  const BrachyPhantomHit& operator = (const BrachyPhantomHit &right);
  int operator == (const BrachyPhantomHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Draw();
  void Print();

private:
  const G4LogicalVolume* logicalVolume;// Logical volume where the hit happened
  G4int xHitPosition; // Hit x coordinate 
  G4int zHitPosition; // Hit z coordinate 
  G4int yHitPosition; // Hit y coordinate 
  G4ThreeVector hitPosition; //position of the hit 
  G4RotationMatrix rotation; // rotation of the logical volume
  G4double energyDeposit; // energy deposit associated with the hit

public:
  //...
  // Set Hit position
  inline void SetCellID(G4int XID,G4int YID,G4int ZID)
  {xHitPosition = XID; zHitPosition = ZID;  yHitPosition = YID;  }
  
  inline G4int GetXID() // Get hit x coordinate 
  {return xHitPosition;}

  inline G4int GetZID() // Get hit z coordinate   
  {return zHitPosition;}

  inline G4int GetYID() // Get hit y coordinate  
  {return yHitPosition;}
   
  inline void SetEdep(G4double edep) //Set hit energy deposit
  {energyDeposit = edep;}

  inline void AddEdep(G4double edep) // Add energy deposit
  {energyDeposit += edep;}

  inline G4double GetEdep() // Get energy deposit
  {return energyDeposit;}

  inline void SetPos(G4ThreeVector xyz) // Set hit position
  {hitPosition = xyz;}

  inline G4ThreeVector GetPos()  // Get hit position
  {return hitPosition;}

  inline void SetRot(G4RotationMatrix rmat) //set rotation
  {rotation = rmat;}

  inline G4RotationMatrix GetRot() //get rotation
  {return rotation;}

  //...
  // It returns the logical volume where the hit happened
  inline const G4LogicalVolume * GetLogicalVolume()
  {return logicalVolume;}
};

typedef G4THitsCollection<BrachyPhantomHit> BrachyPhantomHitsCollection;
extern G4Allocator<BrachyPhantomHit> BrachyPhantomHitAllocator;

inline void* BrachyPhantomHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) BrachyPhantomHitAllocator.MallocSingle();
  return aHit;
}

inline void BrachyPhantomHit::operator delete(void *aHit)
{
  BrachyPhantomHitAllocator.FreeSingle((BrachyPhantomHit*) aHit);
}
#endif


