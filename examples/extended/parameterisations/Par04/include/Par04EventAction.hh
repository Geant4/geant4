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
#ifndef PAR04EVENTACTION_HH
#define PAR04EVENTACTION_HH

#include <G4Types.hh>            // for G4int, G4double
#include <vector>                // for vector
#include "G4Timer.hh"            // for G4Timer
#include "G4UserEventAction.hh"  // for G4UserEventAction
class G4Event;
class Par04DetectorConstruction;

/**
 * @brief Event action class for hits' analysis.
 *
 * Analysis of single-particle events and developed showers in the detector.
 * At the end of the event basic variables are calculated and saved in the
 * histograms.
 * Additionally ntuple with cell energies and IDs (in cylindrical coordinates) is stored.
 *
 */

class Par04EventAction : public G4UserEventAction
{
 public:
  Par04EventAction(Par04DetectorConstruction* aDetector);
  virtual ~Par04EventAction();

  /// Timer is started
  virtual void BeginOfEventAction(const G4Event* aEvent) final;
  /// Hits collection is retrieved, analysed, and histograms are filled.
  virtual void EndOfEventAction(const G4Event* aEvent) final;
  inline std::vector<G4double>& GetCalEdep() { return fCalEdep; }
  inline std::vector<G4int>& GetCalRho() { return fCalRho; }
  inline std::vector<G4int>& GetCalPhi() { return fCalPhi; }
  inline std::vector<G4int>& GetCalZ() { return fCalZ; }

 private:
  /// ID of a hit collection to analyse
  G4int fHitCollectionID;
  /// Timer measurement
  G4Timer fTimer;
  /// Pointer to detector construction to retrieve (once) the detector
  /// dimensions and size of readout
  Par04DetectorConstruction* fDetector = nullptr;
  /// Size of cell along Z axis
  G4double fCellSizeZ = 0;
  /// Size of cell along radius of cylinder
  G4double fCellSizeRho = 0;
  /// Size of cell in azimuthal angle
  G4double fCellSizePhi = 0;
  /// Number of readout cells along radius
  G4int fCellNbRho = 0;
  /// Number of readout cells in azimuthal angle
  G4int fCellNbPhi = 0;
  /// Number of readout cells along z axis
  G4int fCellNbZ = 0;
  /// Cell energy deposits to be stored in ntuple
  std::vector<G4double> fCalEdep;
  /// Cell ID of radius to be stored in ntuple
  std::vector<G4int> fCalRho;
  /// Cell ID of azimuthal angle to be stored in ntuple
  std::vector<G4int> fCalPhi;
  /// Cell ID of z axis to be stored in ntuple
  std::vector<G4int> fCalZ;
};

#endif /* PAR04EVENTACTION_HH */
