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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// $Id: Tst52TrackerHit.hh,v 1.1.2.1 2007-12-10 16:33:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst52TrackerHit_h
#define Tst52TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Tst52TrackerHit : public G4VHit
{
  public:

      Tst52TrackerHit();
     ~Tst52TrackerHit();
      Tst52TrackerHit(const Tst52TrackerHit&);
      const Tst52TrackerHit& operator=(const Tst52TrackerHit&);
      G4int operator==(const Tst52TrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  // Store the energy deposit and position in the hit
      void SetEdep    (G4double de)      { edep = de; };
      void SetPos     (G4int voxel){ voxel_hit = voxel; };
     
      G4double GetEdep()    { return edep; };      
      G4int GetPos(){ return voxel_hit; };
      
  private:
      G4double      edep;
      G4int voxel_hit;
};

typedef G4THitsCollection<Tst52TrackerHit> Tst52TrackerHitsCollection;

extern G4Allocator<Tst52TrackerHit> Tst52TrackerHitAllocator;

inline void* Tst52TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Tst52TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void Tst52TrackerHit::operator delete(void *aHit)
{
  Tst52TrackerHitAllocator.FreeSingle((Tst52TrackerHit*) aHit);
}

#endif
