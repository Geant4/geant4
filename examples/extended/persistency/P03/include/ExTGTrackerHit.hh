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
//
// $Id: ExTGTrackerHit.hh,v 1.2 2010-11-05 08:52:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------------

#ifndef ExTGTrackerHit_h
#define ExTGTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class ExTGTrackerHit : public G4VHit
{
  public:

    ExTGTrackerHit();
   ~ExTGTrackerHit();
    ExTGTrackerHit(const ExTGTrackerHit&);
    const ExTGTrackerHit& operator=(const ExTGTrackerHit&);
    G4int operator==(const ExTGTrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    void Draw();
    void Print();

    void SetTrackID  (G4int track)      { trackID = track; };
    void SetChamberNb(G4int chamb)      { chamberNb = chamb; };  
    void SetEdep     (G4double de)      { edep = de; };
    void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      
    G4int GetTrackID()    { return trackID; };
    G4int GetChamberNb()  { return chamberNb; };
    G4double GetEdep()    { return edep; };      
    G4ThreeVector GetPos(){ return pos; };
      
  private:
  
      G4int         trackID;
      G4int         chamberNb;
      G4double      edep;
      G4ThreeVector pos;
};

// ---------------------------------------------------------------------------

typedef G4THitsCollection<ExTGTrackerHit> ExTGTrackerHitsCollection;

extern G4Allocator<ExTGTrackerHit> ExTGTrackerHitAllocator;

// ---------------------------------------------------------------------------

inline void* ExTGTrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExTGTrackerHitAllocator.MallocSingle();
  return aHit;
}

// ---------------------------------------------------------------------------

inline void ExTGTrackerHit::operator delete(void *aHit)
{
  ExTGTrackerHitAllocator.FreeSingle((ExTGTrackerHit*) aHit);
}

// ---------------------------------------------------------------------------

#endif
