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
// $Id: PhotInCalorHit.hh,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  G4double GetEdep() const { return Edep; }
  G4double GetTrak() const { return TrackLength; }
  G4int   GetNStep() const { return nSteps; }
    
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


