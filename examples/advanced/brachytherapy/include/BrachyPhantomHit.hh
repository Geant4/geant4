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
// $Id: BrachyPhantomHit.hh,v 1.6 2005/11/22 12:47:35 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
//    ********************************
//    *                              *
//    *    BrachyPhantomHit.hh       *
//    *                              *
//    ********************************
//
//Code developed by: S. Guatelli, M.Tropeano
//
// It manages the hits and the enrgy deposit in the phantom associated with 
// each hit ...

#ifndef BrachyPhantomHit_h
#define BrachyPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class BrachyPhantomHit : public G4VHit
{
public:
  BrachyPhantomHit();
  ~BrachyPhantomHit();
  BrachyPhantomHit(const BrachyPhantomHit &right);
  const BrachyPhantomHit& operator = (const BrachyPhantomHit &right);
  int operator == (const BrachyPhantomHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Draw();
  void Print();

private:
  G4double xHitPosition; // Hit x coordinate 
  G4double zHitPosition; // Hit z coordinate 
  G4double yHitPosition; // Hit y coordinate 
  G4double energyDeposit; // energy deposit associated with the hit

public:
  inline G4double GetXID() // Get hit x coordinate 
  {return xHitPosition;}

  inline G4double GetZID() // Get hit z coordinate   
  {return zHitPosition;}

  inline G4double GetYID() // Get hit y coordinate  
  {return yHitPosition;}
   
  inline void SetEdep(G4double edep) //Set hit energy deposit
  {energyDeposit = edep;}

  inline G4double GetEdep() // Get energy deposit
  {return energyDeposit;}

  inline void SetPosition(G4double xx, G4double yy, G4double zz) 
  // Set hit position
  { xHitPosition = xx; yHitPosition = yy; zHitPosition = zz;}
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


