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
// $Id: TargetHit.hh,v 1.1 2003-05-27 13:44:45 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef TargetHit_h
#define TargetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"


class TargetHit : public G4VHit
{
 public:

   TargetHit();
   ~TargetHit();
   TargetHit(const TargetHit&);
   const TargetHit& operator=(const TargetHit&);
   int operator==(const TargetHit&) const;

   inline void* operator new(size_t);
   inline void  operator delete(void*);

   void Draw();
   void Print();
      
 public:
  
   void StoreTheta(G4double cost) {cosTheta = cost;};
   void StoreEnergy(G4double ke) {Ekin = ke;};
   void StoreCharge(G4double q) {Charge = q;};
                 
   G4double GetTheta() {return cosTheta;};
   G4double GetEkin() {return Ekin;};
   G4double GetCharge() {return Charge;};
    
 private:
  
   G4double cosTheta;
   G4double Ekin;     
   G4double Charge;     
};


typedef G4THitsCollection<TargetHit> TargetHitsCollection;

extern G4Allocator<TargetHit> TargetHitAllocator;


inline void* TargetHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) TargetHitAllocator.MallocSingle();
  return aHit;
}


inline void TargetHit::operator delete(void* aHit)
{
  TargetHitAllocator.FreeSingle((TargetHit*) aHit);
}

#endif


