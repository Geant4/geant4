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
//    *    ThyroidHit.hh             *
//    *                              *
//    ********************************

#ifndef ThyroidHit_h
#define ThyroidHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class ThyroidHit : public G4VHit
{
 public:
      ThyroidHit(G4LogicalVolume* logVol);
      ~ThyroidHit();
      ThyroidHit(const ThyroidHit &right);
      const ThyroidHit& operator=(const ThyroidHit &right);
      int operator==(const ThyroidHit &right) const;

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
      G4int m_YID;
      G4int m_ZID;
 
 public:
      inline void SetCellID(G4int XID,G4int YID,G4int ZID)
 	{m_XID = XID;m_YID;m_ZID = ZID;}
      inline G4int GetXID() 
	{return m_XID;}
      inline G4int GetYID() 
	{return m_YID;}
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

typedef G4THitsCollection<ThyroidHit> ThyroidHitsCollection;
extern G4Allocator<ThyroidHit> ThyroidHitAllocator;

inline void* ThyroidHit::operator new(size_t)
{
 void *aHit;
 aHit = (void *) ThyroidHitAllocator.MallocSingle();
 return aHit;
}

inline void ThyroidHit::operator delete(void *aHit)
{
 ThyroidHitAllocator.FreeSingle((ThyroidHit*) aHit);
}

#endif


