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
// $Id: G4MPImanager.hh,v 1.2 2010-05-18 06:03:27 kmura Exp $
// $Name: not supported by cvs2svn $
//
// ====================================================================
//   G4MPImanager.hh
//
//                                         2007 Q
// ====================================================================
#ifndef G4MPI_MANAGER_H
#define G4MPI_MANAGER_H

#include "mpi.h"
#include "globals.hh"
#include <pthread.h>
#include <fstream>

// ====================================================================
//
// class definition
//
// ====================================================================
class G4MPImessenger;
class G4MPIsession;
class G4MPIstatus;
class G4VMPIseedGenerator;

class G4MPImanager {
private:
  static G4MPImanager* theManager;
  G4MPImessenger* messenger;
  G4MPIsession* session;

  G4int verbose;
  G4MPIstatus* status; // status for each node

  // MPI rank
  G4bool isMaster;
  G4bool isSlave;
  G4int rank;
  G4int size;

  // MPI communicator
  MPI::Intracomm COMM_G4COMMAND;

  // cout/cerr control
  G4bool qfcout;
  std::ofstream fscout;

  // init/macro file
  G4bool qinitmacro;
  G4String initFileName;
  G4bool qbatchmode;
  G4String macroFileName;

  // internal use
  void Initialize();
  void ParseArguments(G4int argc, char** argv);
  void Wait(G4int ausec) const;
  void UpdateStatus();

  // for beamOn
  pthread_t threadID;

  // seed generator
  G4VMPIseedGenerator* seedGenerator;

  // parallel parameters
  G4double masterWeight;

public:

  enum { // MPI master rank
    RANK_MASTER= 0
  };

  enum { // MPI tag
    TAG_G4COMMAND= 100,
    TAG_G4STATUS= 200,
    TAG_G4SEED= 300,
    TAG_DATA= 1000
  };

  G4MPImanager();
  G4MPImanager(int argc, char** argv);
  ~G4MPImanager();

  static G4MPImanager* GetManager();

  // set/get methods
  G4MPIsession* GetMPIsession() const;

  G4int GetVerbose() const;
  void SetVerbose(G4int iverbose);

  G4int GetSize() const;
  G4int GetRank() const;

  G4bool IsMaster() const;
  G4bool IsSlave() const;

  G4bool IsInitMacro() const;
  const G4String& GetInitFileName() const;

  G4bool IsBatchMode() const;
  const G4String& GetMacroFileName() const;

  void SetMasterWeight(G4double aweight);
  G4double GetMasterWeight() const;

  G4VMPIseedGenerator* GetSeedGenerator() const;

  // MPI methods
  G4String BcastCommand(const G4String& command);
  void ShowStatus();
  void ShowSeeds();
  void SetSeed(G4int inode, G4long seed);
  void WaitBeamOn();

  // methods for MPI environment
  void DistributeSeeds();
  void ExecuteMacroFile(const G4String& fname, G4bool qbatch=false);
  G4bool CheckThreadStatus();
  void ExecuteThreadCommand(const G4String& command);
  void ExecuteBeamOnThread(const G4String& command);
  void JoinBeamOnThread();

  void BeamOn(G4int nevent, G4bool qdivide=true);
  void Print(const G4String& message);

  // misc
  void ShowHelp() const;

};

// ====================================================================
//   inline functions
// ====================================================================

inline G4MPIsession* G4MPImanager::GetMPIsession() const
{
  return session;
}

inline G4int G4MPImanager::GetVerbose() const
{
  return verbose;
}

inline void G4MPImanager::SetVerbose(G4int iverbose)
{
  G4int lv=iverbose;
  if(iverbose>1) lv=1;
  if(iverbose<0) lv=0;

  verbose= lv;
  return;
}

inline G4int G4MPImanager::GetRank() const
{
  return rank;
}

inline G4int G4MPImanager::GetSize() const
{
  return size;
}

inline G4bool G4MPImanager::IsMaster() const
{
  return isMaster;
}

inline G4bool G4MPImanager::IsSlave() const
{
  return isSlave;
}

inline G4bool G4MPImanager::IsInitMacro() const
{
  return qinitmacro;

}

inline const G4String& G4MPImanager::GetInitFileName() const
{
  return initFileName;

}

inline G4bool G4MPImanager::IsBatchMode() const
{
  return qbatchmode;
}

inline const G4String& G4MPImanager::GetMacroFileName() const
{
  return macroFileName;
}

inline void G4MPImanager::SetMasterWeight(G4double aweight)
{
  masterWeight= aweight;

  if(aweight<0.) masterWeight= 0.;
  if(aweight>1.) masterWeight= 1.;
}

inline G4double G4MPImanager::GetMasterWeight() const
{
  return masterWeight;
}

inline G4VMPIseedGenerator* G4MPImanager::GetSeedGenerator() const
{
  return seedGenerator;
}

#endif
