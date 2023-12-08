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
#ifndef PAR04ACTIONINITIALISATION_HH
#define PAR04ACTIONINITIALISATION_HH

#include "G4VUserActionInitialization.hh"

class Par04DetectorConstruction;
class Par04ParallelFullWorld;

/**
 * @brief Initialization of user actions.
 *
 * Initialises the primary generator, and user actions (event, run) to perform
 * analysis and store histograms.
 *
 */

class Par04ActionInitialisation : public G4VUserActionInitialization
{
 public:
  Par04ActionInitialisation(Par04DetectorConstruction* aDetector,
                            Par04ParallelFullWorld* aParallel);
  ~Par04ActionInitialisation();
  /// Create all user actions.
  virtual void Build() const final;
  /// Create run action in the master thread to allow analysis merging.
  virtual void BuildForMaster() const final;

 private:
  /// Pointer to detector to be passed to event and run actions in order to
  /// retrieve detector dimensions
  Par04DetectorConstruction* fDetector = nullptr;
  Par04ParallelFullWorld* fParallel = nullptr;
};

#endif /* PAR04ACTIONINITIALISATION_HH */
