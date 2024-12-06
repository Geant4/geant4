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
// G4ScoreQuantityMessenger
//
// Class description:
//
// This is a concrete class of G4UImessenger which handles the
// commands for G4ScoringManager.
//
// Author: Tsukasa Aso (KEK), September 2007
// --------------------------------------------------------------------
#ifndef G4ScoreQuantityMessenger_h
#define G4ScoreQuantityMessenger_h 1

#include "G4String.hh"
#include "G4UImessenger.hh"

#include <vector>

class G4ScoringManager;
class G4VScoringMesh;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcommand;

using G4TokenVec = std::vector<G4String>;

class G4ScoreQuantityMessenger : public G4UImessenger
{
 public:

  G4ScoreQuantityMessenger(G4ScoringManager* SManager);

  ~G4ScoreQuantityMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String newValues) override;

  G4String GetCurrentValue(G4UIcommand*) override;

 protected:

  void FillTokenVec(const G4String& newValues, G4TokenVec& token);

  void FParticleCommand(G4VScoringMesh* mesh, G4TokenVec& token);
  void FParticleWithEnergyCommand(G4VScoringMesh* mesh, G4TokenVec& token);

  G4bool CheckMeshPS(G4VScoringMesh* mesh, const G4String& psname,
                     G4UIcommand* command);

 private:

  void QuantityCommands();
  void FilterCommands();

 private:

  G4ScoringManager* fSMan;
  //
  // Quantity commands
  G4UIdirectory* quantityDir;
  G4UIcmdWithAString* qTouchCmd;
  G4UIcmdWithoutParameter* qGetUnitCmd;
  G4UIcmdWithAString* qSetUnitCmd;
  //
  G4UIcommand* qCellChgCmd;
  G4UIcommand* qCellFluxCmd;
  G4UIcommand* qPassCellFluxCmd;
  G4UIcommand* qeDepCmd;
  G4UIcommand* qdoseDepCmd;
  G4UIcommand* qnOfStepCmd;
  G4UIcommand* qnOfSecondaryCmd;
  //
  G4UIcommand* qTrackLengthCmd;
  G4UIcommand* qPassCellCurrCmd;
  G4UIcommand* qPassTrackLengthCmd;
  G4UIcommand* qFlatSurfCurrCmd;
  G4UIcommand* qFlatSurfFluxCmd;
  G4UIcommand* qVolFluxCmd;
  G4UIcommand* qNofCollisionCmd;
  G4UIcommand* qPopulationCmd;
  G4UIcommand* qTrackCountCmd;
  G4UIcommand* qTerminationCmd;
  G4UIcommand* qMinKinEAtGeneCmd;
  G4UIcommand* qStepCheckerCmd;

  //
  // Filter commands
  G4UIdirectory* filterDir;
  G4UIcmdWithAString* fchargedCmd;
  G4UIcmdWithAString* fneutralCmd;
  G4UIcommand* fkinECmd;
  G4UIcommand* fparticleCmd;
  G4UIcommand* fparticleKinECmd;
  //
};

#endif
