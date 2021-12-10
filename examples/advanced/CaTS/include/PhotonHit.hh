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

// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PhotonHit.hh
/// \brief Definition of the CaTS::PhotonHit class

#pragma once

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class PhotonHit : public G4VHit
{
 public:
  PhotonHit();
  ~PhotonHit() = default;
  PhotonHit(const PhotonHit&);
  const PhotonHit& operator=(const PhotonHit&);
  G4bool operator==(const PhotonHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;

  inline void Print() final
  {
    G4cout << "PhotonHit id: " << fid << " pid: " << fpid
           << " wavelength: " << fwavelength << " time: " << ftime
           << " position X: " << fposition.getX() << " Y: " << fposition.getY()
           << " Z: " << fposition.getZ()
           << " direction X: " << fdirection.getX()
           << " Y: " << fdirection.getY() << " Z: " << fdirection.getZ()
           << " polarization: X:" << fpolarization.getX()
           << " Y: " << fpolarization.getY() << " Z: " << fpolarization.getZ()
           << G4endl;
  }

  PhotonHit(unsigned id, unsigned pid, G4double wavelength, G4double time,
            G4ThreeVector position, G4ThreeVector direction,
            G4ThreeVector polarization);
  inline void SetWavelength(G4double wavelength) { fwavelength = wavelength; }
  inline G4double GetWavelength() { return fwavelength; }
  inline void SetPolarization(G4ThreeVector polarization)
  {
    fpolarization = polarization;
  }
  inline G4ThreeVector GetPolarization() const { return fpolarization; }
  inline void SetDirection(G4ThreeVector direction) { fdirection = direction; }
  inline G4ThreeVector GetDirection() const { return fdirection; }
  inline void SetPosition(G4ThreeVector position) { fposition = position; }
  inline G4ThreeVector GetPosition() const { return fposition; }
  inline void SetTime(G4double time) { ftime = time; }
  inline G4double GetTime() const { return ftime; }
  inline void SetPid(unsigned pid) { fpid = pid; }
  inline unsigned GetPid() const { return fpid; }
  inline void SetId(unsigned id) { fid = id; }
  inline unsigned GetId() const { return fid; }

 private:
  unsigned int fid{ 0 };
  unsigned int fpid{ 0 };
  G4double fwavelength{ 0 };
  G4double ftime{ 0 };
  G4ThreeVector fposition{ 0, 0, 0 };
  G4ThreeVector fdirection{ 0, 0, 0 };
  G4ThreeVector fpolarization{ 0, 0, 0 };
};

using PhotonHitsCollection = G4THitsCollection<PhotonHit>;
extern G4ThreadLocal G4Allocator<PhotonHit>* PhotonHitAllocator;

inline void* PhotonHit::operator new(size_t)
{
  if(!PhotonHitAllocator)
  {
    PhotonHitAllocator = new G4Allocator<PhotonHit>;
  }
  return (void*) PhotonHitAllocator->MallocSingle();
}

inline void PhotonHit::operator delete(void* aHit)
{
  PhotonHitAllocator->FreeSingle((PhotonHit*) aHit);
}
