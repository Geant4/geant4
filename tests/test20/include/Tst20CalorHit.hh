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
// $Id: Tst20CalorHit.hh,v 1.4 2007-11-09 18:32:59 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $


#ifndef Tst20CalorHit_h
#define Tst20CalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"


class Tst20CalorHit : public G4VHit
{
public:

  Tst20CalorHit();
  ~Tst20CalorHit();
  Tst20CalorHit(const Tst20CalorHit&);
  const Tst20CalorHit& operator=(const Tst20CalorHit&);
  bool operator==(const Tst20CalorHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Print();
      
  void AddEnergyDeposit(G4double energy, G4double length);
                 
  G4double GetEnergyDeposit() { return energyDeposit; }
  G4double GetTrackLength() { return trackLength; }
     
private:
  
  G4double energyDeposit;
  G4double trackLength;
      
};


typedef G4THitsCollection<Tst20CalorHit> Tst20CalorHitsCollection;

extern G4Allocator<Tst20CalorHit> Tst20CalorHitAllocator;


inline void* Tst20CalorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) Tst20CalorHitAllocator.MallocSingle();
  return aHit;
}


inline void Tst20CalorHit::operator delete(void* aHit)
{
  Tst20CalorHitAllocator.FreeSingle((Tst20CalorHit*) aHit);
}

#endif


