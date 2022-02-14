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

#ifndef SILICONPIXELHIT_HH
#define SILICONPIXELHIT_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4Types.hh"

#include <vector>

/**
 * @brief Silicon pixel hit
 *
 * Stores information of energy deposited in the silicon pixel.
 *
 * Hits can be digitised, to take into account the time window for the deposits
 * as well as to set the time information based on the energy threshold.
 * By default no time window is used (all energy deposits are counted), as well
 * as no energy threshold is used (time of the first energy deposit, however
 * small, is counted as hit time).
 *
 */

class SiliconPixelHit : public G4VHit {
public:
  ///	Constructor
  /// @param[in] aName Name of the pixel volume
  /// @param[in] aCopyNumSensor ID of the sensor
  /// @param[in] aCopyNumCell ID of the cell
  SiliconPixelHit(G4String aName, G4int aCopyNumSensor, G4int aCopyNumCell);
  ~SiliconPixelHit(){};
  /// Draw pixels
  void Draw();
  /// Get hit ID calculated as 1000 * sensorID + cellID
  G4int ID() { return 1000 * fCopyNumSensor + fCopyNumCell; }
  /// Add non-zero energy deposit to vector of deposits
  /// @param[in] aEnergy Deposited energy
  /// @param[in] aTime Time of deposit
  inline void AddEdep(const G4double aEnergy, const G4double aTime) {
    if (aEnergy > 0)
      fEdep.push_back(std::make_pair(aEnergy, aTime));
  }
  /// Add non-zero non-ionizing energy deposit to vector of deposits
  /// @param[in] aEnergy Deposited energy
  /// @param[in] aTime Time of deposit
  inline void AddEdepNonIonizing(const G4double aEnergy, const G4double aTime) {
    if (aEnergy > 0)
      fEdepNonIonizing.push_back(std::make_pair(aEnergy, aTime));
  }
  /// Digitise hit
  /// Calculate time of hit as global time of energy deposit which added
  /// to hit energy exceeds the energy threshold. Take into account only
  /// deposits with global time within the time window.
  /// @param[in] aTimeWindow Maximal global time for deposit, caounted from
  /// the time of the first deposit
  /// @param[in] aToaThreshold Energy threshold, first deposit that adds to
  /// the hit energy and exceeds the threshold is counted as time of
  /// arrival.
  void Digitise(const G4double aTimeWindow, const G4double aToaThreshold);
  /// Set hit position
  /// @param[in] x X position
  /// @param[in] y Y position
  /// @param[in] z Z position
  inline void SetPosition(G4double x, G4double y, G4double z) {
    fPosX = x;
    fPosY = y;
    fPosZ = z;
  }
  /// Get hit X position
  inline G4double GetX() const { return fPosX; }
  /// Get hit Y position
  inline G4double GetY() const { return fPosY; }
  /// Get hit Z position
  inline G4double GetZ() const { return fPosZ; }
  /// Check if hit is valid
  inline G4bool isValidHit() const { return fIsValidHit; }
  /// Get hit energy
  inline G4double GetEdep() const { return fEdepDigi; }
  /// Get hit non-ionizing energy
  inline G4double GetEdepNonIonizing() const { return fEdepNonIonizingDigi; }
  /// Get time of arrival
  inline G4double GetTOA() const { return fTimeOfArrival; }
  /// Get time of arrival from the last energy deposit
  inline G4double GetLastTOA() const { return fTimeOfArrivalLast; }

private:
  ///	Name of the logical volume
  G4String fVolumeName = "";
  /// ID of the sensor
  G4int fCopyNumCell = -1;
  /// ID of the cell
  G4int fCopyNumSensor = -1;
  /// Position along x axis
  G4double fPosX = -1;
  /// Position along y axis
  G4double fPosY = -1;
  /// Position along z axis
  G4double fPosZ = -1;
  /// Vector of energy deposits (and their global time)
  std::vector<std::pair<G4double, G4double>> fEdep;
  /// Vector of non-ionizing energy deposits (and their global time)
  std::vector<std::pair<G4double, G4double>> fEdepNonIonizing;
  /// Flag indicating if hit is valid (digitised and with non-zero energy)
  G4bool fIsValidHit = false;
  /// Energy of the digitised hit
  G4double fEdepDigi = -1;
  /// Non-ionizing energy of the digitised hit
  G4double fEdepNonIonizingDigi = -1;
  /// Time of arrival of the digitised hit
  G4double fTimeOfArrival = -1;
  /// Last time of arrival of the digitised hit
  G4double fTimeOfArrivalLast = -1;
};

typedef G4THitsCollection<SiliconPixelHit> SiliconPixelHitCollection;

#endif /* SILICONPIXELHIT_HH */