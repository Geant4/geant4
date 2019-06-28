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
/// \file persistency/P01/include/ExP01TrackerHit.hh
/// \brief Definition of the ExP01TrackerHit class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExP01TrackerHit_h
#define ExP01TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Hit implementation for the persistency example

class ExP01TrackerHit : public G4VHit
{
  public:

      ExP01TrackerHit();
     ~ExP01TrackerHit();
      ExP01TrackerHit(const ExP01TrackerHit&);
      const ExP01TrackerHit& operator=(const ExP01TrackerHit&);
      G4bool operator==(const ExP01TrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      virtual void Draw();
      virtual void Print();

  public:
  
      void SetTrackID  (G4int track)      { fTrackID = track; };
      void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };  
      void SetEdep     (G4double de)      { fEdep = de; };
      void SetPos      (G4ThreeVector xyz){ fPos = xyz; };
      
      G4int GetTrackID()    { return fTrackID; };
      G4int GetChamberNb()  { return fChamberNb; };
      G4double GetEdep()    { return fEdep; };      
      G4ThreeVector GetPos(){ return fPos; };
      
  private:
  
      G4int         fTrackID;
      G4int         fChamberNb;
      G4double      fEdep;
      G4ThreeVector fPos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using ExP01TrackerHitsCollection = G4THitsCollection<ExP01TrackerHit>;

extern G4ThreadLocal G4Allocator<ExP01TrackerHit>* ExP01TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* ExP01TrackerHit::operator new(size_t)
{
  if (!ExP01TrackerHitAllocator)
  {
    ExP01TrackerHitAllocator = new G4Allocator<ExP01TrackerHit>;
  }
  return (void *) ExP01TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void ExP01TrackerHit::operator delete(void *aHit)
{
  ExP01TrackerHitAllocator->FreeSingle((ExP01TrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
