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
/// @file G4MPIstatus.cc
/// @brief status of MPI application

#include "G4ApplicationState.hh"
#include "G4MPIstatus.hh"

// --------------------------------------------------------------------------
G4MPIstatus::G4MPIstatus()
  : rank_(0), run_id_(0), nevent_to_be_processed_(0), event_id_(0),
    cputime_(0.), g4state_(G4State_Quit)
{
  timer_ = new G4Timer;
}

// --------------------------------------------------------------------------
G4MPIstatus::~G4MPIstatus()
{
  delete timer_;
}

// --------------------------------------------------------------------------
void G4MPIstatus::StartTimer()
{
  timer_-> Start();
}

// --------------------------------------------------------------------------
void G4MPIstatus::StopTimer()
{
  timer_-> Stop();
}

// --------------------------------------------------------------------------
void G4MPIstatus::SetStatus(G4int arank, G4int runid, G4int noe, G4int evtid,
                            G4ApplicationState state)
{
  rank_ = arank;
  run_id_ = runid;
  nevent_to_be_processed_ = noe;
  event_id_ = evtid;
  g4state_ = state;
  if ( timer_-> IsValid() ) cputime_= timer_-> GetRealElapsed();
  else cputime_ = 0.;
}

// --------------------------------------------------------------------------
void G4MPIstatus::Pack(G4int* data) const
{
  data[0] = rank_;
  data[1] = run_id_;
  data[2] = nevent_to_be_processed_;
  data[3] = event_id_;
  data[4] = g4state_;

  G4double* ddata = (G4double*)(data+5);
  ddata[0] = cputime_;
}

// --------------------------------------------------------------------------
void G4MPIstatus::UnPack(G4int* data)
{
  rank_ = data[0];
  run_id_ = data[1];
  nevent_to_be_processed_ = data[2];
  event_id_ = data[3];
  g4state_ = (G4ApplicationState)data[4];

  G4double* ddata = (G4double*)(data+5);
  cputime_ = ddata[0];
}

// --------------------------------------------------------------------------
void G4MPIstatus::Print() const
{
  // * rank= 001 run= 10002 event= 00001 / 100000 state= Idle"
  G4cout << "* rank= " << rank_
         << " run= " << run_id_
         << " event= " << event_id_ << " / " << nevent_to_be_processed_
         << " state= " << GetStateString(g4state_)
         << " time= " << cputime_ << "s"
         << G4endl;
}

// --------------------------------------------------------------------------
G4String G4MPIstatus::GetStateString(G4ApplicationState astate) const
{
  G4String sname;

  switch(astate) {
  case G4State_PreInit:
    sname = "PreInit";
    break;
  case G4State_Init:
    sname = "Init";
    break;
  case G4State_Idle:
    sname = "Idle";
    break;
  case G4State_GeomClosed:
    sname = "GeomClosed";
    break;
  case G4State_EventProc:
    sname = "EventProc";
    break;
  case G4State_Quit:
    sname = "Quit";
    break;
  case G4State_Abort:
    sname = "Abort";
    break;
  default:
    sname = "Unknown";
    break;
  }

  return sname;
}
