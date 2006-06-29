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
// $Id: RE01TrackerHit.hh,v 1.3 2006-06-29 17:43:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


#ifndef RE01TrackerHit_h
#define RE01TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class RE01TrackerHit : public G4VHit
{
  public:

      RE01TrackerHit();
      ~RE01TrackerHit();
      RE01TrackerHit(const RE01TrackerHit &right);
      const RE01TrackerHit& operator=(const RE01TrackerHit &right);
      G4int operator==(const RE01TrackerHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;
      G4int trackID;

  public:
      inline void SetEdep(G4double de)
      { edep = de; }
      inline G4double GetEdep() const
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos() const
      { return pos; }
      inline void SetTrackID(G4int i)
      { trackID = i; }
      inline G4int GetTrackID() const
      { return trackID; }
};

typedef G4THitsCollection<RE01TrackerHit> RE01TrackerHitsCollection;

extern G4Allocator<RE01TrackerHit> RE01TrackerHitAllocator;

inline void* RE01TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) RE01TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void RE01TrackerHit::operator delete(void *aHit)
{
  RE01TrackerHitAllocator.FreeSingle((RE01TrackerHit*) aHit);
}

#endif
