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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//
// Code based on the basic example B02

#ifndef SensitiveDetectorHit_h
#define SensitiveDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh" // FOR MT

class SensitiveDetectorHit : public G4VHit
{
  public:
    SensitiveDetectorHit();
    SensitiveDetectorHit(const SensitiveDetectorHit&);
    virtual ~SensitiveDetectorHit();

    // operators
    const SensitiveDetectorHit& operator=(const SensitiveDetectorHit&);
    G4bool operator==(const SensitiveDetectorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetEdep     (G4double de)      { fEdep = de; };
   
    G4double GetEdep() const     { return fEdep; };
 
  private:

      G4double      fEdep;
};

typedef G4THitsCollection<SensitiveDetectorHit> SensitiveDetectorHitsCollection;

extern G4ThreadLocal G4Allocator<SensitiveDetectorHit>* SensitiveDetectorHitAllocator;

inline void* SensitiveDetectorHit::operator new(size_t)
{
  if(!SensitiveDetectorHitAllocator)
      SensitiveDetectorHitAllocator = new G4Allocator<SensitiveDetectorHit>;
  return (void *) SensitiveDetectorHitAllocator->MallocSingle();
}

inline void SensitiveDetectorHit::operator delete(void *hit)
{
  SensitiveDetectorHitAllocator->FreeSingle((SensitiveDetectorHit*) hit);
}

#endif
