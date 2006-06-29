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
// $Id: Tst14CalorHit.hh,v 1.7 2006-06-29 21:40:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef Tst14CalorHit_h
#define Tst14CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"


class Tst14CalorHit : public G4VHit
{
  public:

      Tst14CalorHit();
     ~Tst14CalorHit();
      Tst14CalorHit(const Tst14CalorHit&);
      const Tst14CalorHit& operator=(const Tst14CalorHit&);
      int operator==(const Tst14CalorHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Print();
      
  public:
  
      void AddAbs(G4double de, G4double dl) {EdepAbs += de; TrackLengthAbs += dl;};
      void AddGap(G4double de, G4double dl) {EdepGap += de; TrackLengthGap += dl;};      
                 
      G4double GetEdepAbs()     { return EdepAbs; };
      G4double GetTrakAbs()     { return TrackLengthAbs; };
      G4double GetEdepGap()     { return EdepGap; };
      G4double GetTrakGap()     { return TrackLengthGap; };
     
  private:
  
      G4double EdepAbs, TrackLengthAbs;
      G4double EdepGap, TrackLengthGap;
      
};

typedef G4THitsCollection<Tst14CalorHit> Tst14CalorHitsCollection;

extern G4Allocator<Tst14CalorHit> Tst14CalorHitAllocator;


inline void* Tst14CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) Tst14CalorHitAllocator.MallocSingle();
  return aHit;
}


inline void Tst14CalorHit::operator delete(void* aHit)
{
  Tst14CalorHitAllocator.FreeSingle((Tst14CalorHit*) aHit);
}

#endif


