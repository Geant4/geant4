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
// $Id: ExN07CalorHit.hh,v 1.1 2003/03/10 01:43:35 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifndef ExN07CalorHit_h
#define ExN07CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class ExN07CalorHit : public G4VHit
{
 public:
   ExN07CalorHit();
   virtual ~ExN07CalorHit();
   ExN07CalorHit(const ExN07CalorHit&);
   const ExN07CalorHit& operator=(const ExN07CalorHit&);
   int operator==(const ExN07CalorHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   virtual void Draw();
   virtual void Print();
      
 public:
   inline void AddEnergy(G4double de)
   { Edep += de; }
   inline void AddStep(G4double dl)
   {
     TrackLength += dl; 
     nStep++;
   }
   inline G4double GetEdep() const
   { return Edep; }
   inline G4double GetTrak() const
   { return TrackLength; }
   inline G4int GetNStep() const
   { return nStep; }
    
 private:
   G4double Edep;
   G4double TrackLength;
   G4int    nStep;
};

typedef G4THitsCollection<ExN07CalorHit> ExN07CalorHitsCollection;

extern G4Allocator<ExN07CalorHit> ExN07CalorHitAllocator;

inline void* ExN07CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) ExN07CalorHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN07CalorHit::operator delete(void* aHit)
{
  ExN07CalorHitAllocator.FreeSingle((ExN07CalorHit*) aHit);
}


#endif


