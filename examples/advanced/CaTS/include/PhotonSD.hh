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
/// \file PhotonSD.hh
/// \brief Definition of the CaTS::PhotonSD class

#pragma once

#include "G4VSensitiveDetector.hh"
#include "G4PhysicalConstants.hh"
#include "PhotonHit.hh"
#include <G4String.hh>
#include <G4Types.hh>
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PhotonSD : public G4VSensitiveDetector
{
 public:
  PhotonSD(G4String);
  ~PhotonSD() = default;
  void Initialize(G4HCofThisEvent*) final;
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) final;
#ifdef WITH_G4OPTICKS
  void AddOpticksHits();
#endif
  void EndOfEvent(G4HCofThisEvent* hitCollection) final;

 private:
  PhotonHitsCollection* fPhotonHitsCollection{ 0 };
  G4int fHCID{ 0 };
  static constexpr G4double hc = (h_Planck * c_light) / (CLHEP::eV * CLHEP::nm);
  G4bool verbose{ false };
  inline G4double etolambda(G4double E)
  {
    // input photon energy in eV
    // return wavelength in nm:
    // lambda = h c/e
    //
    return hc / E;
  }
};
