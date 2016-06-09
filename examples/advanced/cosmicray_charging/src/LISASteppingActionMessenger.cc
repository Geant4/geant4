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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISASteppingActionMessenger class                                *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************

#include "LISASteppingActionMessenger.hh"

#include <sstream>

#include "LISASteppingAction.hh"


LISASteppingActionMessenger::LISASteppingActionMessenger
  (LISASteppingAction* stepAct) : steppingAction(stepAct){


  newDir = new G4UIdirectory("/surveys/");
  newDir->SetGuidance("Particle survey control commands.");

  // set flag for spectrum saving
  SetFlagSpectrum = new G4UIcmdWithABool("/surveys/surveyTestMasses",this);
  SetFlagSpectrum->SetGuidance("Dump energies (eV) for e- and hadrons");
  SetFlagSpectrum->SetGuidance("entering and leaving the test masses.");
  SetFlagSpectrum->SetGuidance("'electrons_in.out' : e- entering TM");
  SetFlagSpectrum->SetGuidance("'electrons_out.out': e- leaving TM");
  SetFlagSpectrum->SetGuidance("'hadrons_in.out'   : hadrons entering TM");
  SetFlagSpectrum->SetGuidance("'hadrons_out.out'  : hadrons leaving TM");
  SetFlagSpectrum->SetGuidance("Default = false");
  SetFlagSpectrum->SetParameterName("FlagSpectrum", false);
  SetFlagSpectrum->AvailableForStates(G4State_Idle);

}


LISASteppingActionMessenger::~LISASteppingActionMessenger() {

  delete SetFlagSpectrum;  
  delete newDir;

}


void LISASteppingActionMessenger::SetNewValue
   (G4UIcommand* command, G4String newValue) { 

  if(command == SetFlagSpectrum) {
    G4int vl;
    const char* t = newValue;
    std::istringstream is(t);
    is >> vl;
    steppingAction->SetFlagSpectrum(vl!=0);
  }



}

