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
// Rich advanced example for Geant4
// RichTbHit.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbHit_h
#define RichTbHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

class RichTbVisManager;
class RichTbHit : public G4VHit
{
  public:

      RichTbHit();
      virtual ~RichTbHit();
      RichTbHit(const RichTbHit &right);
      const RichTbHit& operator=(const RichTbHit &right);
      int operator==(const RichTbHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

  void Draw();
  void DrawWithVisM(RichTbVisManager*);
  void Print();

  private:
      G4double edep;
      G4ThreeVector posAtSilicon;
      G4ThreeVector posAtPhotoCathode;
      G4int CurHpdNum;
      G4int CurSectNum;
      G4int CurPixelNum;
  public:
  inline void SetEdep(G4double de)      { edep = de; };
  inline G4double GetEdep()      { return edep; };
  inline void SetPos(G4ThreeVector xyz)      { posAtSilicon = xyz; };
  inline G4ThreeVector GetPos()      { return posAtSilicon; };
  inline void SetPosPC(G4ThreeVector xyzPC)  { posAtPhotoCathode = xyzPC; };
  inline G4ThreeVector GetPosPC()      { return posAtPhotoCathode; };
  inline void SetCurHpdNum (G4int ihp ) { CurHpdNum = ihp; } ;
  inline G4int GetCurHpdNum()           { return CurHpdNum ; };
  inline void SetCurSectNum (G4int isector ) { CurSectNum = isector; } ;
  inline G4int GetCurSectNum()           { return CurSectNum ; };
  inline void SetCurPixNum (G4int ipx ) { CurPixelNum = ipx; };
  inline G4int GetCurPixNum()           { return CurPixelNum; };
  inline void AddEdep( G4double addenergy ) { edep += addenergy; }  
};

typedef G4THitsCollection<RichTbHit> RichTbHitsCollection;

extern G4Allocator<RichTbHit> RichTbHitAllocator;

inline void* RichTbHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) RichTbHitAllocator.MallocSingle();
  return aHit;
}

inline void RichTbHit::operator delete(void *aHit)
{
  RichTbHitAllocator.FreeSingle((RichTbHit*) aHit);
}



#endif




