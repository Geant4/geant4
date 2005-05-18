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
// It manages the hits and the enrgy deposit in the phantom associated with 
// each hit ...

#ifndef HadrontherapyPhantomHit_h
#define HadrontherapyPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class HadrontherapyPhantomHit : public G4VHit
{
public:
  HadrontherapyPhantomHit();
  ~HadrontherapyPhantomHit();
  HadrontherapyPhantomHit(const HadrontherapyPhantomHit &right);
  const HadrontherapyPhantomHit& operator = (const HadrontherapyPhantomHit &right);
  int operator == (const HadrontherapyPhantomHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Draw();
  void Print();

private:
  G4int xHitID; // Hit x voxel 
  G4int zHitID; // Hit z voxel
  G4int yHitID; // Hit y voxel 
  G4double energyDeposit; // energy deposit associated with the hit

public:
  //...
  // Set Hit position
  inline void SetCellID(G4int XID,G4int YID,G4int ZID)
  {xHitID = XID; zHitID = ZID;  yHitID = YID;  }
  
  inline G4int GetXID() // Get hit x coordinate 
  {return xHitID;}

  inline G4int GetZID() // Get hit z coordinate   
  {return zHitID;}

  inline G4int GetYID() // Get hit y coordinate  
  {return yHitID;}
   
  inline void SetEdep(G4double edep) //Set hit energy deposit
  {energyDeposit = edep;}

  inline void SetEdepAndPosition(G4int xx, G4int yy, G4int zz, G4double eDep)
  {
    xHitID = xx;
    yHitID = yy;
    zHitID = zz;
    energyDeposit = eDep;
  }

  inline void AddEdep(G4double edep) // Add energy deposit
  {energyDeposit += edep;}

  inline G4double GetEdep() // Get energy deposit
  {return energyDeposit;}
};

typedef G4THitsCollection<HadrontherapyPhantomHit> HadrontherapyPhantomHitsCollection;
extern G4Allocator<HadrontherapyPhantomHit> HadrontherapyPhantomHitAllocator;

inline void* HadrontherapyPhantomHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HadrontherapyPhantomHitAllocator.MallocSingle();
  return aHit;
}

inline void HadrontherapyPhantomHit::operator delete(void *aHit)
{
  HadrontherapyPhantomHitAllocator.FreeSingle((HadrontherapyPhantomHit*) aHit);
}
#endif


