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
/// @file G4MPIstatus.hh
/// @brief status of MPI application

#ifndef G4MPI_STATUS_H
#define G4MPI_STATUS_H

#include "globals.hh"
#include "G4ApplicationState.hh"
#include "G4Timer.hh"

class G4Timer;

class G4MPIstatus {
public:
  G4MPIstatus();
  ~G4MPIstatus();

  // set/get functions
  void SetStatus(G4int arank, G4int runid, G4int noe, G4int evtid,
                 G4ApplicationState state);

  G4int GetRank() const;
  G4int GetRunID() const;
  G4int GetNEventToBeProcessed() const;
  G4int GetEventID() const;
  G4double GetCPUTime() const;
  G4ApplicationState GetG4State() const;

  // for pickling
  G4int SizeOf() const;
  void Pack(G4int* data) const;
  void UnPack(G4int* data);

  // for timer
  void StartTimer();
  void StopTimer();

  void Print() const;

  enum { kNSIZE = 10 };

private:
  G4int rank_;
  G4int run_id_;
  G4int nevent_to_be_processed_;
  G4int event_id_;
  G4double cputime_;
  G4ApplicationState g4state_;
  G4Timer* timer_;

  G4String GetStateString(G4ApplicationState astate) const;
};

// ====================================================================
inline G4int G4MPIstatus::SizeOf() const
{
  return kNSIZE;
}

inline G4int G4MPIstatus::GetRank() const
{
  return rank_;
}

inline G4int G4MPIstatus::GetRunID() const
{
  return run_id_;
}

inline G4int G4MPIstatus::GetNEventToBeProcessed() const
{
  return nevent_to_be_processed_;
}

inline G4int G4MPIstatus::GetEventID() const
{
  return event_id_;
}

inline G4double G4MPIstatus::GetCPUTime() const
{
  return cputime_;
}

inline G4ApplicationState G4MPIstatus::GetG4State() const
{
  return g4state_;
}

#endif
