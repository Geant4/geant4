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

// Author: Ivana Hrivnacova, 18/10/2022  (ivana@ipno.in2p3.fr)

#include "G4GenericAnalysisMessenger.hh"
#include "G4GenericAnalysisManager.hh"

#include "G4UIcmdWithAString.hh"

//_____________________________________________________________________________
G4GenericAnalysisMessenger::G4GenericAnalysisMessenger(G4GenericAnalysisManager* manager)
  : fManager(manager)
{
  fSetDefaultFileTypeCmd = CreateCommand<G4UIcmdWithAString>(
    "setDefaultFileType", "Set default output file type", "DefaultFileType", false);
#ifdef TOOLS_USE_HDF5
  fSetDefaultFileTypeCmd->SetCandidates("csv hdf5 root xml");
#else
  fSetDefaultFileTypeCmd->SetCandidates("csv root xml");
#endif
}

//_____________________________________________________________________________
G4GenericAnalysisMessenger::~G4GenericAnalysisMessenger() = default;

//
//
// public functions
//

//_____________________________________________________________________________
void G4GenericAnalysisMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetDefaultFileTypeCmd.get() ) {
    fManager->SetDefaultFileType(newValues);
    return;
  }
}
