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
///  \file XSHistoManagerMessenger.hh
///  \brief UI commands for XS study.
//
//  Author: G.Hugo, 06 January 2023
//
// ***************************************************************************
//
//      XSHistoManagerMessenger
//
///  UI commands for XS study.
//
//  NB: Note that it is not necessary to have POINTERS to the UI commands.
//
// ***************************************************************************

#ifndef HISTO_MANAGER_MESSENGER_HH
#define HISTO_MANAGER_MESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


class XSHistoManager;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class XSHistoManagerMessenger: public G4UImessenger {

public:

  XSHistoManagerMessenger(XSHistoManager* const histoManager);

  virtual void SetNewValue(G4UIcommand*, G4String) override;


private:

  XSHistoManager* fHisto = nullptr;

  G4UIcmdWithAString fOutputFileNameCmd;
  G4UIcmdWithAString fParticleNameCmd;
  G4UIcmdWithAString fElementNameCmd;
  G4UIcmdWithAString fNonElementaryMaterialNameCmd;
  G4UIcmdWithAnInteger fNumBinsCmd;
  G4UIcmdWithADoubleAndUnit fMinKineticEnergyCmd;
  G4UIcmdWithADoubleAndUnit fMaxKineticEnergyCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
