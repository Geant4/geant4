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
// Class G4PSVolumeFlux
//
// Class description:
//  Scorer that scores number of tracks coming into the associated volume.
//  Optionally number of tracks can be divided by the surface area to
//  score the volume current and divided by cos(theta) for volume flux
//  where theta is the incident angle
//
//  - Created   M. Asai, Sept. 2020

#ifndef G4PSVolumeFlux_h
#define G4PSVolumeFlux_h 1

#include "G4VPrimitivePlotter.hh"
#include "G4THitsMap.hh"
#include "G4PSDirectionFlag.hh"

class G4PSVolumeFlux : public G4VPrimitivePlotter
{
 public:
  G4PSVolumeFlux(G4String name, G4int direction = 1, G4int depth = 0);
  ~G4PSVolumeFlux() override = default;

 public:
  void Initialize(G4HCofThisEvent*) override;
  void clear() override;
  void PrintAll() override;

  void SetDivAre(G4bool val) { divare = val; }
  void SetDivCos(G4bool val) { divcos = val; }
 
 protected:
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

 private:
  G4int HCID{-1};
  G4int fDirection;
  G4THitsMap<G4double>* EvtMap{nullptr};
  G4bool divare = false;  // divide by the surface area
  G4bool divcos = false;  // divide by cos(theta)
};

#endif
