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
/// \file TimeStepAction.hh
/// \brief Definition of the TimeStepAction class

#ifndef ITACTION_H
#define ITACTION_H

#include "G4MolecularConfiguration.hh"
#include "G4MoleculeIterator.hh"
#include "G4MoleculeTable.hh"
#include "G4UserTimeStepAction.hh"

class DetectorConstruction;
class G4Molecule;

class TimeStepAction final : public G4UserTimeStepAction {
public:
  TimeStepAction();
  ~TimeStepAction() override = default;
  TimeStepAction(const TimeStepAction &other) = delete;
  TimeStepAction &operator=(const TimeStepAction &other) = delete;
  void UserPreTimeStepAction() override;
  void EndProcessing() override { ; }
  static void Save(const G4MolecularConfiguration *molconf);
  static void SaveMoleculeInfo(const G4Track *track, G4int molID, const G4String &
                               /*moleculeName*/);
private:
  const DetectorConstruction *fpDetector = nullptr;
};

#endif  // ITACTION_H