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
//
/// \file FlashStackingAction.hh
/// \brief Definition of the FlashStackingAction class

#ifndef FlashStackingAction_h
#define FlashStackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include <fstream>

class FlashStackingAction : public G4UserStackingAction {
public:
  FlashStackingAction();
  virtual ~FlashStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

private:
  std::ofstream OpticFile;
  G4int fScintillationCounter;
  G4int fCerenkovCounter;
  G4int fBremstralung;
  G4int fFluo;
  G4double Annihil_x;
  G4double Annihil_y;
  G4double Annihil_z;
};

#endif
