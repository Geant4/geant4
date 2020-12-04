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
#ifndef SILICONPIXELSD_HH
#define SILICONPIXELSD_HH

#include "G4VSensitiveDetector.hh"
#include "SiliconPixelHit.hh"
#include "G4Types.hh"

#include <map>

/**
 * @brief Sensitive detector for silicon pixels
 *
 * Processes step information and stores it in silicon hits.
 * Information is stored in units of:
 * - keV for energy
 * - cm for position
 *
 */

class SiliconPixelSD : public G4VSensitiveDetector {
public:
  explicit SiliconPixelSD(G4String name);
  ~SiliconPixelSD();
  /// Hits are processed and added to the temporary map of hit ID to hit
  /// pointer. Energy is stored in units of keV, and position in cm.
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *ROhist);

  void Initialize(G4HCofThisEvent *HCE);
  /// Temporary map of hits is stored in hit collection, to be retrieved
  /// for analysis by the event action
  void EndOfEvent(G4HCofThisEvent *HCE);

private:
  /// Hit collection stored in the event, filled in at the end of event based
  /// on temporary hits
  SiliconPixelHitCollection *fHitCollection = nullptr;
  /// ID of hit collection
  G4int fHCID = -1;
  /// Temporary map of hits (ID: hit) collected within one event
  std::map<G4int, SiliconPixelHit *> fTmpHits;
};
#endif /* SILICONPIXELSD_HH */