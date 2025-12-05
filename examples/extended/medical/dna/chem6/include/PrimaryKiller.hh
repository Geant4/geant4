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
/// \file PrimaryKiller.hh
/// \brief Definition of the PrimaryKiller class

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#ifndef CHEM6_PrimaryKiller_h
#define CHEM6_PrimaryKiller_h 1

#include <G4THitsMap.hh>
#include <G4UImessenger.hh>
#include <G4VPrimitiveScorer.hh>
#include <G4UnitsTable.hh>

class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;

/** \file PrimaryKiller.hh*/

// Description:
//   Kill the primary particle:
//   - either after a given energy loss
//   - or after the primary particle has reached a given energy

class PrimaryKiller final : public G4VPrimitiveScorer, public G4UImessenger {
public:
  explicit PrimaryKiller(const G4String &name, G4int depth = 0);

  ~PrimaryKiller() override = default;

  /** Set energy under which the particle should be
   killed*/
  inline void SetEnergyThreshold(const G4double energy) { fKineticE_Min = energy; }

  /** Set the energy loss from which the primary is
   killed*/
  inline void SetMinLossEnergyLimit(double energy) { fELossRange_Min = energy; }

  /** Set the energy loss from which the event is
   aborted*/
  inline void SetMaxLossEnergyLimit(double energy) { fELossRange_Max = energy; }

  /** Method related to G4UImessenger
      used to control energy cuts through macro file
   */
  void SetNewValue(G4UIcommand *command, G4String newValue) override;

  void Initialize(G4HCofThisEvent *) override {fELoss = 0.;};

  void EndOfEvent(G4HCofThisEvent *) override {};


protected:
  G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;

private:
  G4double fELoss = 0; // cumulated energy loss by the primary
  G4double fELossRange_Min = DBL_MAX; // fELoss from which the primary is killed
  G4double fELossRange_Max = DBL_MAX; // fELoss from which the event is aborted
  G4double fKineticE_Min = 0; // kinetic energy below which the primary is killed
  G4ThreeVector fPhantomSize = {1 * CLHEP::km, 1 * CLHEP::km, 1 * CLHEP::km};

  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpELossUI;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpAbortEventIfELossUpperThan;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fpMinKineticE;
  std::unique_ptr<G4UIcmdWith3VectorAndUnit> fpSizeUI;
};
#endif
