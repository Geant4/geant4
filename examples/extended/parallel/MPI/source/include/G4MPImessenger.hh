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

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class G4MPImanager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcommand;

class G4MPImessenger : public G4UImessenger {
public:
  G4MPImessenger();
  ~G4MPImessenger();

  virtual void SetNewValue(G4UIcommand* command, G4String newValue);
  virtual G4String GetCurrentValue(G4UIcommand* command);

  void SetTargetObject(G4MPImanager* mpi_manager);

private:
  DISALLOW_COPY_AND_ASSIGN(G4MPImessenger);

  G4MPImanager* g4mpi_;

  // /mpi
  G4UIdirectory* dir_;

  G4UIcmdWithAnInteger* verbose_;
  G4UIcmdWithoutParameter* status_;

  G4UIcmdWithAString* execute_;

  G4UIcommand* beam_on_;
  G4UIcommand* dot_beam_on_;
  G4UIcmdWithADouble* master_weight_;

  G4UIcmdWithoutParameter* show_seeds_;
  G4UIcmdWithAnInteger* set_master_seed_;
  G4UIcommand* set_seed_;
};

// ====================================================================
inline void G4MPImessenger::SetTargetObject(G4MPImanager* mpi_manager)
{
  g4mpi_ = mpi_manager;
}

#endif
