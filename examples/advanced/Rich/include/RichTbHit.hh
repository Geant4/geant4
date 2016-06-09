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

class G4VVisManager;
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
  void DrawWithVisM(G4VVisManager*);
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




