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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef HadrontherapyDetectorHit_h
#define HadrontherapyDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class HadrontherapyDetectorHit : public G4VHit
{
public:
  HadrontherapyDetectorHit();
  ~HadrontherapyDetectorHit();
 
  HadrontherapyDetectorHit(const HadrontherapyDetectorHit &right);
 
  const HadrontherapyDetectorHit& operator = (const HadrontherapyDetectorHit &right);
 
  int operator == (const HadrontherapyDetectorHit &right) const;

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

typedef G4THitsCollection<HadrontherapyDetectorHit> HadrontherapyDetectorHitsCollection;
extern G4Allocator<HadrontherapyDetectorHit> HadrontherapyDetectorHitAllocator;

inline void* HadrontherapyDetectorHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HadrontherapyDetectorHitAllocator.MallocSingle();
  return aHit;
}

inline void HadrontherapyDetectorHit::operator delete(void *aHit)
{
  HadrontherapyDetectorHitAllocator.FreeSingle((HadrontherapyDetectorHit*) aHit);
}
#endif


