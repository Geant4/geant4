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
// $Id: PhotInCalorHit.hh,v 1.4 2006/06/29 16:24:35 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#ifndef PhotInCalorHit_h
#define PhotInCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class PhotInCalorHit : public G4VHit
{
public:
  PhotInCalorHit();
  virtual ~PhotInCalorHit();
  PhotInCalorHit(const PhotInCalorHit&);
  const PhotInCalorHit& operator=(const PhotInCalorHit&);
  int operator==(const PhotInCalorHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  virtual void Draw();
  virtual void Print();
      
  void AddEnergy(G4double de) { Edep += de; }

  void AddStep(G4double dl)
  {
    TrackLength += dl; 
    nSteps++;
  }
  void ResetEDepos() { Edep=0.; }
  void ResetTrackL() { TrackLength=0.; }
  void ResetNSteps() { nSteps=0; }

  G4double GetEDepos() const { return Edep; }
  G4double GetTrackL() const { return TrackLength; }
  G4int    GetNSteps() const { return nSteps; }
    
private: //--- BODY ---
  G4double Edep;
  G4double TrackLength;
  G4int    nSteps;
};

typedef G4THitsCollection<PhotInCalorHit> PhotInCalorHitsCollection;

extern G4Allocator<PhotInCalorHit> PhotInCalorHitAllocator;

inline void* PhotInCalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) PhotInCalorHitAllocator.MallocSingle();
  return aHit;
}

inline void PhotInCalorHit::operator delete(void* aHit)
  { PhotInCalorHitAllocator.FreeSingle((PhotInCalorHit*) aHit); }


#endif


