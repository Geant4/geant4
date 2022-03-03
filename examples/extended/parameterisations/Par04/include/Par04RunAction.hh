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
#ifndef PAR04RUNACTION_HH
#define PAR04RUNACTION_HH

#include "G4UserRunAction.hh"
#include "Par04PrimaryGeneratorAction.hh"
#include <CLHEP/Units/SystemOfUnits.h>       // for GeV
#include <G4String.hh>                       // for G4String
#include <G4ThreeVector.hh>                  // for G4ThreeVector
#include <G4Types.hh>                        // for G4int
#include <G4VUserPrimaryGeneratorAction.hh>  // for G4VUserPrimaryGeneratorA...
#include <string>                            // for basic_string
#include "G4Event.hh"                        // for G4Event
#include "G4ParticleGun.hh"                  // for G4ParticleGun
#include "G4ParticleTable.hh"                // for G4ParticleTable
#include "G4SystemOfUnits.hh"                // for GeV
#include "Par04EventInformation.hh"          // for Par04EventInformation
class G4ParticleDefinition;
class Par04EventAction;
class G4Run;
class Par04DetectorConstruction;

/**
 * @brief Run action
 *
 * Create analysis file and define control histograms for showers in detectors.
 * Histograms are configured taking into account the dimensions of the readout mesh.
 * Ntuple with hits is also stored. It contains energy of the primary particle,
 * coordinates (cylindrical) of the hit and the deposited energy.
 *
 */

class Par04RunAction : public G4UserRunAction
{
 public:
  /// Constructor. Defines the histograms.
  Par04RunAction(Par04DetectorConstruction* aDetector, Par04EventAction* aEventAction);
  virtual ~Par04RunAction();

  /// Open the file for the analysis
  virtual void BeginOfRunAction(const G4Run*) final;
  /// Write and close the file
  virtual void EndOfRunAction(const G4Run*) final;

 private:
  /// Pointer to detector construction to retrieve the detector dimensions to
  /// setup the histograms
  Par04DetectorConstruction* fDetector;
  /// Pointer to event action to save hits
  Par04EventAction* fEventAction;
};

#endif /* PAR04RUNACTION_HH */
