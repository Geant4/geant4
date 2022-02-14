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
#ifndef PAR03EVENTACTION_HH
#define PAR03EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Timer.hh"

class Par03DetectorConstruction;

/**
 * @brief Event action class for hits' analysis.
 *
 * Analysis of single-particle events and developed showers in the detector.
 * At the end of the event basic variables are calculated and saved in the
 * histograms.
 *
 */

class Par03EventAction : public G4UserEventAction
{
 public:
  Par03EventAction(Par03DetectorConstruction* aDetector);
  virtual ~Par03EventAction();

  /// Timer is started
  virtual void BeginOfEventAction(const G4Event* aEvent) final;
  /// Hits collection is retrieved, analysed, and histograms are filled.
  virtual void EndOfEventAction(const G4Event* aEvent) final;

 private:
  /// ID of a hit collection to analyse
  G4int fHitCollectionID;
  /// Timer measurement
  G4Timer fTimer;
  /// Pointer to detector construction to retrieve (once) the detector
  /// dimensions
  Par03DetectorConstruction* fDetector;
  /// Size of cell along Z axis
  G4double fCellSizeZ = 0;
  /// Size of cell along radius of cylinder
  G4double fCellSizeRho = 0;
};

#endif /* PAR03EVENTACTION_HH */
