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
//    ********************************
//    *                              *
//    *    BrachyWaterBoxHit.hh      *
//    *                              *
//    ********************************

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
  BrachyPhantomHit(G4LogicalVolume* logVol,G4int XID,G4int ZID);
  ~BrachyPhantomHit();
  BrachyPhantomHit(const BrachyPhantomHit &right);
  const BrachyPhantomHit& operator=(const BrachyPhantomHit &right);
  int operator==(const BrachyPhantomHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Draw();
  void Print();

private:
  G4ThreeVector m_Pos;
  G4RotationMatrix m_Rot;
  const G4LogicalVolume* m_pLogV;

  G4double m_Edep;

  G4int m_XID;
  G4int m_ZID;
 
public:
  inline void SetCellID(G4int XID,G4int ZID)
  {m_XID = XID;m_ZID = ZID;}
  inline G4int GetXID() 
  {return m_XID;}
  inline G4int GetZID() 
  {return m_ZID;}
  inline void SetEdep(G4double edep)
  {m_Edep = edep;}
  inline void AddEdep(G4double edep)
  {m_Edep += edep;}
  inline G4double GetEdep()
  {return m_Edep;}
  inline void SetPos(G4ThreeVector xyz)
  {m_Pos = xyz;}
  inline G4ThreeVector GetPos()
  {return m_Pos;}
  inline void SetRot(G4RotationMatrix rmat)
  {m_Rot = rmat;}
  inline G4RotationMatrix GetRot()
  {return m_Rot;}
  inline const G4LogicalVolume * GetLogV()
  {return m_pLogV;}
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


