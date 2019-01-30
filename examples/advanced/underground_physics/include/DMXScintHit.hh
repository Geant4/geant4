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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef DMXScintHit_h
#define DMXScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DMXScintHit : public G4VHit
{
  public:

      DMXScintHit();
      ~DMXScintHit();
      DMXScintHit(const DMXScintHit&);
      const DMXScintHit& operator=(const DMXScintHit&);
      G4bool operator==(const DMXScintHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetEdep           (G4double de)       { edep = de; };
      void SetPos            (G4ThreeVector xyz) { pos = xyz; };
      void SetParticle       (G4String name)     { particleName = name; };
      void SetParticleEnergy (G4double e1)       { particleEnergy = e1; };
      void SetTime           (G4double t2)       { time = t2; };


      G4double GetEdep()                         { return edep; };      
      G4ThreeVector GetPos()                     { return pos; };
      G4String GetParticle()                     { return particleName;};
      G4double GetParticleEnergy()               { return particleEnergy;};
      G4double GetTime()                         { return time; };      


  private:
      G4double      edep;
      G4ThreeVector pos;
      G4double      time;
      G4String      particleName;
      G4double      particleEnergy;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<DMXScintHit> DMXScintHitsCollection;

extern G4ThreadLocal G4Allocator<DMXScintHit> *DMXScintHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* DMXScintHit::operator new(size_t)
{
  if (!DMXScintHitAllocator)
    DMXScintHitAllocator = new G4Allocator<DMXScintHit>;
  return (void*) DMXScintHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void DMXScintHit::operator delete(void* aHit)
{
  DMXScintHitAllocator->FreeSingle((DMXScintHit*) aHit);
}

#endif

