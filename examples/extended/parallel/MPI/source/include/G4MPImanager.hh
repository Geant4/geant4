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
/// @file G4MPImanager.hh
/// @brief MPI manager class

#ifndef G4MPI_MANAGER_H
#define G4MPI_MANAGER_H

#include "mpi.h"
#include <fstream>
#include <pthread.h>
#include "globals.hh"

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class G4MPImessenger;
class G4MPIsession;
class G4MPIstatus;
class G4VMPIseedGenerator;
class G4VMPIextraWorker;

class G4MPImanager {
public:
  // MPI master rank
  enum { kRANK_MASTER = 0 };

  enum { // MPI tag
    kTAG_G4COMMAND = 100,
    kTAG_G4STATUS = 200,
    kTAG_G4SEED = 300,
    kTAG_DATA = 1000,
    kTAG_HISTO = 1001,
    kTAG_RUN = 1002,
    kTAG_CMDSCR = 1003,
    kTAG_NTUPLE = 1004
  };

  G4MPImanager(int nof_extra_workers = 0);
  G4MPImanager(int argc, char** argv, int nof_extra_workers = 0);
  ~G4MPImanager();

  static G4MPImanager* GetManager();

  // set/get methods
  G4MPIsession* GetMPIsession() const;

  G4int GetVerbose() const;
  void SetVerbose(G4int iverbose);

  G4int GetTotalSize() const;    // get size of all ranks
  G4int GetActiveSize() const;   // get size of ranks wher RunBeamOn is called
  G4int GetRank() const;

  G4bool IsMaster() const;
  G4bool IsSlave() const;
  G4bool IsExtraWorker() const;

  G4bool IsInitMacro() const;
  const G4String& GetInitFileName() const;

  G4bool IsBatchMode() const;
  const G4String& GetMacroFileName() const;

  void SetMasterWeight(G4double aweight);
  G4double GetMasterWeight() const;

  void SetExtraWorker(G4VMPIextraWorker* extraWorker);
  G4VMPIextraWorker* GetExtraWorker() const;

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

  const MPI::Intracomm* GetComm() const { return &COMM_G4COMMAND_; }
  const MPI_Comm* GetProcessingComm() const { return &processing_comm_; }
  const MPI_Comm* GetCollectingComm() const { return &collecting_comm_; }
  const MPI_Comm* GetAllComm() const { return &all_comm_; }
private:
  DISALLOW_COPY_AND_ASSIGN(G4MPImanager);

  // internal use
  void Initialize();
  void ParseArguments(G4int argc, char** argv);
  void UpdateStatus();

  static G4MPImanager* g4mpi_;
  G4MPImessenger* messenger_;
  G4MPIsession* session_;
  G4VMPIextraWorker* extra_worker_;

  // seed generator
  G4VMPIseedGenerator* seed_generator_;

  G4MPIstatus* status_; // status for each node

  G4int verbose_;

  // MPI rank
  G4bool is_master_;
  G4bool is_slave_;
  G4bool is_extra_worker_;
  G4int rank_;
  G4int size_;  // processing comm size
  G4int world_size_;  // world comm size

  // MPI communicator (when no extra ranks)
  MPI::Intracomm COMM_G4COMMAND_;
  // MPI communicator (processing ranks - if ntuple merging)
  MPI_Comm processing_comm_;
  // MPI communicator (collecting ranks - if ntuple merging)
  MPI_Comm collecting_comm_;
  // MPI communicator (all ranks - if ntuple mergins)
  MPI_Comm all_comm_;
  // Interim data - need to be freed
  MPI_Group world_group_;
  MPI_Group processing_group_;
  MPI_Group collecting_group_;
  MPI_Group all_group_;

  // cout/cerr control
  G4bool qfcout_;
  std::ofstream fscout_;

  // init/macro file
  G4bool qinitmacro_;
  G4String init_file_name_;
  G4bool qbatchmode_;
  G4String macro_file_name_;

  // for beamOn
  pthread_t thread_id_;

  // parallel parameters
  G4double master_weight_;
  G4int nof_extra_workers_;
};

// ====================================================================
inline G4MPIsession* G4MPImanager::GetMPIsession() const
{
  return session_;
}

inline G4int G4MPImanager::GetVerbose() const
{
  return verbose_;
}

inline void G4MPImanager::SetVerbose(G4int iverbose)
{
  G4int lv = iverbose;
  if( iverbose > 1 ) lv = 1;
  if( iverbose < 0 ) lv = 0;

  verbose_ = lv;
  return;
}

inline G4int G4MPImanager::GetRank() const
{
  return rank_;
}

inline G4int G4MPImanager::GetTotalSize() const
{
  return world_size_;
}

inline G4int G4MPImanager::GetActiveSize() const
{
  return size_;
}

inline G4bool G4MPImanager::IsMaster() const
{
  return is_master_;
}

inline G4bool G4MPImanager::IsSlave() const
{
  return is_slave_;
}

inline G4bool G4MPImanager::IsExtraWorker() const
{
  return is_extra_worker_;
}

inline G4bool G4MPImanager::IsInitMacro() const
{
  return qinitmacro_;

}

inline const G4String& G4MPImanager::GetInitFileName() const
{
  return init_file_name_;

}

inline G4bool G4MPImanager::IsBatchMode() const
{
  return qbatchmode_;
}

inline const G4String& G4MPImanager::GetMacroFileName() const
{
  return macro_file_name_;
}

inline void G4MPImanager::SetMasterWeight(G4double aweight)
{
  master_weight_ = aweight;

  if( aweight < 0. ) master_weight_ = 0.;
  if( aweight > 1. ) master_weight_ = 1.;
}

inline G4double G4MPImanager::GetMasterWeight() const
{
  return master_weight_;
}

inline G4VMPIextraWorker* G4MPImanager::GetExtraWorker() const
{
  return extra_worker_;
}

inline G4VMPIseedGenerator* G4MPImanager::GetSeedGenerator() const
{
  return seed_generator_;
}

#endif
