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
/// @file G4MPImessenger.hh
/// @brief Define MPI commands

#ifndef G4MPI_MESSENGER_H
#define G4MPI_MESSENGER_H

#include "G4UImessenger.hh"

class G4MPImanager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcommand;

class G4MPImessenger : public G4UImessenger {
public:
  G4MPImessenger(G4MPImanager* manager);
  ~G4MPImessenger();

  virtual void SetNewValue(G4UIcommand* command, G4String newValue);
  virtual G4String GetCurrentValue(G4UIcommand* command);

private:
  G4MPImanager* g4MPI;

  // /mpi
  G4UIdirectory* dir;

  G4UIcmdWithAnInteger* verbose;
  G4UIcmdWithoutParameter* status;

  G4UIcmdWithAString* execute;

  G4UIcommand* beamOn;
  G4UIcommand* dotbeamOn;
  G4UIcmdWithADouble* masterWeight;

  G4UIcmdWithoutParameter* showSeeds;
  G4UIcmdWithAnInteger* setMasterSeed;
  G4UIcommand* setSeed;
};

#endif
