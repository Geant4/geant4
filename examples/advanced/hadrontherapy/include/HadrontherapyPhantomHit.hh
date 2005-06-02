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
// $Id: HadrontherapyPhantomHit.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyPhantomHit_h
#define HadrontherapyPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

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

private:
  G4int xHitID; // Hit x voxel 
  G4int zHitID; // Hit z voxel
  G4int yHitID; // Hit y voxel 
  G4double energyDeposit; // Energy deposit associated with the hit

public:
  // Methods to get the information - energy deposit and associated
  // position in the phantom - of the hits stored in the hits collection  
 
  inline G4int GetXID() // Get x index of the voxel 
  {return xHitID;}

  inline G4int GetZID() // Get y index of the voxel   
  {return zHitID;}

  inline G4int GetYID() // Get z index of the voxel  
  {return yHitID;}
   
  inline G4double GetEdep() // Get energy deposit
  {return energyDeposit;}
 
  // Methods to store the information of the hit ( energy deposit, position in the phantom )
  // in the hits collection

  inline void SetEdepAndPosition(G4int xx, G4int yy, G4int zz, G4double eDep)
  {
    xHitID = xx;
    yHitID = yy;
    zHitID = zz;
    energyDeposit = eDep;
  }
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


