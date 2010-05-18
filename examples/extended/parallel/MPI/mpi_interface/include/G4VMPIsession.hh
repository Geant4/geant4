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
// $Id: G4VMPIsession.hh,v 1.2 2010-05-18 06:05:05 kmura Exp $
// $Name: not supported by cvs2svn $
//
// ====================================================================
//   G4VMPIsession.hh
//
//                                         2007 Q
// ====================================================================
#ifndef G4VMPI_SESSION_H
#define G4VMPI_SESSION_H

#include "G4VBasicShell.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

typedef void* (*Func_t)(void *);  // for thead function

class G4MPImanager;
class G4VUIshell;

class G4VMPIsession : public G4VBasicShell {
protected:
  // MPI handling
  G4MPImanager* g4MPI;

  // MPI node info (cache)
  G4bool isMaster;
  G4bool isSlave;
  G4int rank;

  G4int ExecCommand(G4String acommand);
  G4String TruncateCommand(const G4String& command) const;
  G4String BypassCommand(const G4String& command) const;

  // for help operation
  virtual G4bool GetHelpChoice(G4int& aval);
  virtual void ExitHelp();
  
public:
  G4VMPIsession();
  ~G4VMPIsession();

  // concrete implementations
  virtual void PauseSessionStart(G4String msg);

  virtual G4int ReceiveG4cout(G4String coutString);
  virtual G4int ReceiveG4cerr(G4String cerrString);

};

#endif
