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
/// @file G4VMPIsession.hh
/// @brief A base class for MPI sessions

#ifndef G4VMPI_SESSION_H
#define G4VMPI_SESSION_H

#include "G4VBasicShell.hh"

typedef void* (*Func_t)(void *);  // for thread function

class G4MPImanager;
class G4VUIshell;

class G4VMPIsession : public G4VBasicShell {
public:
  G4VMPIsession();
  ~G4VMPIsession();

  virtual void PauseSessionStart(const G4String& msg);

  virtual G4int ReceiveG4cout(const G4String& coutString);
  virtual G4int ReceiveG4cerr(const G4String& cerrString);

protected:
  G4MPImanager* g4mpi_;

  // MPI node info (cache)
  G4bool is_master_;
  G4bool is_slave_;
  G4int rank_;

  G4int ExecCommand(const G4String& acommand);
  G4String TruncateCommand(const G4String& command) const;
  G4String BypassCommand(const G4String& command) const;

  // for help operation
  virtual G4bool GetHelpChoice(G4int& aval);
  virtual void ExitHelp() const;
};

#endif
