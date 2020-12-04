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
#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Types.hh"

#include <vector>

class G4GenericMessenger;

/**
 * @brief Event action class
 *
 * Fills ntuples with information about hits in sensitive detectors.
 *
 */

class EventAction : public G4UserEventAction {
public:
  EventAction();
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event *event);
  virtual void EndOfEventAction(const G4Event *event);

  /// Vector of primary particles in an event: PDG type
  std::vector<G4int> fPrimariesPDG;
  /// Vector of primary particles in an event: particle energy (in GeV)
  std::vector<G4double> fPrimariesEnergy;
  /// Vector of primary particles in an event: vertex position x (in cm)
  std::vector<G4double> fPrimariesX;
  /// Vector of primary particles in an event: vertex position y (in cm)
  std::vector<G4double> fPrimariesY;
  /// Vector of primary particles in an event: vertex position z (in cm)
  std::vector<G4double> fPrimariesZ;

  /// Vector of hits in silicon sensors: hit ID
  std::vector<G4int> fSiHitsID;
  /// Vector of hits in silicon sensors: hit position x (in cm)
  std::vector<G4double> fSiHitsX;
  /// Vector of hits in silicon sensors: hit position y (in cm)
  std::vector<G4double> fSiHitsY;
  /// Vector of hits in silicon sensors: hit position z (in cm)
  std::vector<G4double> fSiHitsZ;
  /// Vector of hits in silicon sensors: hit energy (in keV)
  std::vector<G4double> fSiHitsEdep;
  /// Vector of hits in silicon sensors: hit non-ionizing energy (in keV)
  std::vector<G4double> fSiHitsEdepNonIonising;
  /// Vector of hits in silicon sensors: hit time of arrival (in ns)
  /// calculated as global time of energy deposit which added to hit energy
  /// exceeds the toa threshold (by default threshold is 0, so it is the first
  /// deposit)
  std::vector<G4double> fSiHitsTOA;
  /// Vector of hits in silicon sensors: hit time of last arrival (in ns)
  /// calculated as global time of energy deposit which is the last deposit
  /// that fits within the digitisation time window (by default window is
  /// undefined, so it is the last deposit)
  std::vector<G4double> fSiHitsTOAlast;
  /// Vector of hits in silicon sensors: hit type
  /// Simulation defines only hits of type 0 (hexagonal cells, no calibration or
  /// edge cell is constructed)
  std::vector<G4int> fSiHitsType;

  /// Vector of hits in SiPM: hit ID
  std::vector<G4int> fSiPMhitsID;
  /// Vector of hits in SiPM: hit position x (in cm)
  std::vector<G4double> fSiPMhitsX;
  /// Vector of hits in SiPM: hit position y (in cm)
  std::vector<G4double> fSiPMhitsY;
  /// Vector of hits in SiPM: hit position z (in cm)
  std::vector<G4double> fSiPMhitsZ;
  /// Vector of hits in SiPM: hit energy (in keV)
  std::vector<G4double> fSiPMhitsEdep;
  /// Vector of hits in SiPM: hit non-ionizing energy (in keV)
  std::vector<G4double> fSiPMhitsEdepNonIonising;
  /// Vector of hits in SiPM: hit time of last arrival (in ns)
  /// calculated as global time of energy deposit which is the last deposit
  /// that fits within the digitisation time window (by default window is
  /// undefined, so it is the last deposit)
  std::vector<G4double> fSiPMhitsTOA;
  /// Vector of hits in silicon sensors: hit type
  /// Simulation defines only hits of type 1
  std::vector<G4int> fSiPMhitsType;

private:
  /// Define UI commands: digitisation of hits with the time cut on deposits
  /// (time window), and the energy threshold for the first deposit counted as
  /// the time of arrival
  void DefineCommands();
  /// Pointer to the messenger for UI commands
  G4GenericMessenger *fMessenger = nullptr;
  /// Time window for hit digitisation (in ns)
  /// By default undefined window indicates the last created deposit will set up
  /// the time
  /// Can be changed by the UI command /HGCalTestbeam/hits/timeCut
  G4double fHitTimeCut = -1;
  /// Time of arrival threshold (in keV)
  /// Default value of 0 indicates the first created deposit will set up the
  /// time, independent on the amount of deposited energy
  /// Can be changed by the UI command /HGCalTestbeam/hits/toaThreshold
  G4double fToaThreshold = 0;
};

#endif /* EVENTACTION_HH */
