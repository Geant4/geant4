//    ********************************
//    *                              *
//    *    BrachyWaterBoxHit.hh      *
//    *                              *
//    ********************************

#ifndef BrachyWaterBoxHit_h
#define BrachyWaterBoxHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class BrachyWaterBoxHit : public G4VHit
{
 public:
      BrachyWaterBoxHit(G4LogicalVolume* logVol,G4int XID,G4int ZID);
      ~BrachyWaterBoxHit();
      BrachyWaterBoxHit(const BrachyWaterBoxHit &right);
      const BrachyWaterBoxHit& operator=(const BrachyWaterBoxHit &right);
      int operator==(const BrachyWaterBoxHit &right) const;

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

typedef G4THitsCollection<BrachyWaterBoxHit> BrachyWaterBoxHitsCollection;
extern G4Allocator<BrachyWaterBoxHit> BrachyWaterBoxHitAllocator;

inline void* BrachyWaterBoxHit::operator new(size_t)
{
 void *aHit;
 aHit = (void *) BrachyWaterBoxHitAllocator.MallocSingle();
 return aHit;
}

inline void BrachyWaterBoxHit::operator delete(void *aHit)
{
 BrachyWaterBoxHitAllocator.FreeSingle((BrachyWaterBoxHit*) aHit);
}

#endif


