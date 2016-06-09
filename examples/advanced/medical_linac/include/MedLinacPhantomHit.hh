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
// $Id: MedLinacPhantomHit.hh,v 1.3 2006/06/29 16:03:49 gunter Exp $
//
//
// Code developed by: M. Piergentili
// It manages the hits and the enrgy deposit in the phantom associated with 
// each hit ...

#ifndef MedLinacPhantomHit_h
#define MedLinacPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class MedLinacPhantomHit : public G4VHit
{
public:
  MedLinacPhantomHit(G4LogicalVolume* ,G4int ,G4int ,G4int );
  ~MedLinacPhantomHit();
  MedLinacPhantomHit(const MedLinacPhantomHit &right);
  const MedLinacPhantomHit& operator = (const MedLinacPhantomHit &right);
  int operator == (const MedLinacPhantomHit &right) const;

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

typedef G4THitsCollection<MedLinacPhantomHit> MedLinacPhantomHitsCollection;
extern G4Allocator<MedLinacPhantomHit> MedLinacPhantomHitAllocator;

inline void* MedLinacPhantomHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) MedLinacPhantomHitAllocator.MallocSingle();
  return aHit;
}

inline void MedLinacPhantomHit::operator delete(void *aHit)
{
  MedLinacPhantomHitAllocator.FreeSingle((MedLinacPhantomHit*) aHit);
}
#endif


