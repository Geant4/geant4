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
/// \file ExN04TrackerHit.hh
/// \brief Definition of the ExN04TrackerHit class
//
//
// --------------------------------------------------------------

#ifndef ExN04TrackerHit_h
#define ExN04TrackerHit_h 1


//MSH_include_begin
#include "MarshaledG4String.h"
//MSH_include_end

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//MSH_BEGIN
class ExN04TrackerHit : public G4VHit 
{
  public:

      ExN04TrackerHit();
     ~ExN04TrackerHit();

      ExN04TrackerHit(const ExN04TrackerHit &right);
      const ExN04TrackerHit& operator=(const ExN04TrackerHit &right);
      G4bool operator==(const ExN04TrackerHit &right) const;

      inline void* operator new(size_t);
      inline void  operator delete(void* aHit);

      void Draw();
      void Print();


  private:
  G4ThreeVector pos; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetPos(); }]
  [elementSet: { $THIS->SetPos($ELEMENT); }] */

  G4double edep; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetEdep(); }]
  [elementSet: { $THIS->SetEdep($ELEMENT); }] */


  public:

    inline void SetEdep(G4double de)            {edep=de;}
    inline void SetPos(G4ThreeVector xyz)       {pos=xyz;}

    inline G4double GetEdep()                   const {return edep;}
    inline G4ThreeVector GetPos()               const {return pos;}

};
//MSH_END

// vector collection of one type of hits
typedef G4THitsCollection<ExN04TrackerHit> ExN04TrackerHitsCollection;

extern G4Allocator<ExN04TrackerHit> ExN04TrackerHitAllocator;

inline void* ExN04TrackerHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) ExN04TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04TrackerHit::operator delete(void* aHit)
{
  ExN04TrackerHitAllocator.FreeSingle((ExN04TrackerHit*) aHit);
}

#endif
