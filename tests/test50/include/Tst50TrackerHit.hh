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
// $Id: Tst50TrackerHit.hh,v 1.4 2006-06-29 22:05:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG   first implementation    

#ifndef Tst50TrackerHit_h
#define Tst50TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Tst50TrackerHit : public G4VHit
{
public:

  Tst50TrackerHit();
  ~Tst50TrackerHit();
  Tst50TrackerHit(const Tst50TrackerHit&);
  const Tst50TrackerHit& operator=(const Tst50TrackerHit&);
  int operator==(const Tst50TrackerHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

public:
  void SetEdep     (G4double de)      { edep = de; };
  void AddEnergy    (G4double de)      { edep += de; };    
  G4double GetEdep()    { return edep; };      

private:
  G4double      edep;
};

typedef G4THitsCollection<Tst50TrackerHit> Tst50TrackerHitsCollection;

extern G4Allocator<Tst50TrackerHit> Tst50TrackerHitAllocator;

inline void* Tst50TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Tst50TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void Tst50TrackerHit::operator delete(void *aHit)
{
  Tst50TrackerHitAllocator.FreeSingle((Tst50TrackerHit*) aHit);
}

#endif


