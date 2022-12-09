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
// EmStandardPhysicsTrackingManager
//
// Class description:
//
// An implementation of the G4VTrackingManager interface for e-/e+ and gamma
// with the same processes as G4EmStandardPhysics.
//
// Original author: Jonas Hahnfeld, 2021

#ifndef EmStandardPhysicsTrackingManager_h
#define EmStandardPhysicsTrackingManager_h 1

#include "G4VTrackingManager.hh"
#include "globals.hh"

class G4eMultipleScattering;
class G4CoulombScattering;
class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4ComptonScattering;
class G4GammaConversion;
class G4PhotoElectricEffect;
class G4RayleighScattering;

class EmStandardPhysicsTrackingManager : public G4VTrackingManager
{
 public:
  EmStandardPhysicsTrackingManager();
  ~EmStandardPhysicsTrackingManager();

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  void HandOverOneTrack(G4Track* aTrack) override;

 private:
  void TrackElectron(G4Track* aTrack);
  void TrackPositron(G4Track* aTrack);
  void TrackGamma(G4Track* aTrack);

  struct
  {
    G4eMultipleScattering* msc;
    G4eIonisation* ioni;
    G4eBremsstrahlung* brems;
    G4CoulombScattering* ss;
  } fElectronProcs;

  struct
  {
    G4eMultipleScattering* msc;
    G4eIonisation* ioni;
    G4eBremsstrahlung* brems;
    G4eplusAnnihilation* annihilation;
    G4CoulombScattering* ss;
  } fPositronProcs;

  struct
  {
    G4PhotoElectricEffect* pe;
    G4ComptonScattering* compton;
    G4GammaConversion* conversion;
    G4RayleighScattering* rayleigh;
  } fGammaProcs;

  static EmStandardPhysicsTrackingManager* fMasterTrackingManager;
};

#endif
