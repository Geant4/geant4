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
// $Id: HadrontherapyHit.hh,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------

#ifndef HadrontherapyHit_h
#define HadrontherapyHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

// ------------------------------------------------------------------
class HadrontherapyHit : public G4VHit
{
public:
  HadrontherapyHit(G4int);
  ~HadrontherapyHit();
  HadrontherapyHit(const HadrontherapyHit &right);
  const HadrontherapyHit& operator=(const HadrontherapyHit &right);
  int operator==(const HadrontherapyHit &right) const;
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();
  
private:
     
  G4int slice;
  G4double m_Edep;

public: 
  inline G4int GetSliceID() {return slice;}
  inline void SetSliceID(G4int sliceID)
  {slice= sliceID;}
  inline void SetEdep(G4double edep)
  {m_Edep = edep;}
  inline void AddEdep(G4double edep)
  {m_Edep += edep;}
  inline G4double GetEdep()
  {return m_Edep;}
};

typedef G4THitsCollection<HadrontherapyHit> HadrontherapyHitsCollection;
extern G4Allocator<HadrontherapyHit> HadrontherapyHitAllocator;
inline void* HadrontherapyHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HadrontherapyHitAllocator.MallocSingle();
  return aHit;
}

inline void HadrontherapyHit::operator delete(void *aHit)
{
  HadrontherapyHitAllocator.FreeSingle((HadrontherapyHit*) aHit);
}
#endif


