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
#include "G4VUserMPIrunMerger.hh"

#include "G4MPImanager.hh"
#include "G4MPIutils.hh"

#include <algorithm>
#include <assert.h>
#include <functional>
#include <mpi.h>

G4VUserMPIrunMerger::G4VUserMPIrunMerger(const G4Run* aRun, G4int destination, G4int ver)
  : outputBuffer(nullptr),
    outputBufferSize(0),
    outputBufferPosition(0),
    ownsBuffer(false),
    destinationRank(destination),
    run(const_cast<G4Run*>(aRun)),
    commSize(0),
    verbose(ver),
    bytesSent(0)
{}

#define DMSG(LVL, MSG)         \
  {                            \
    if (verbose > LVL) {       \
      G4cout << MSG << G4endl; \
    }                          \
  }

void G4VUserMPIrunMerger::Send(const unsigned int destination)
{
  assert(run != nullptr);
  G4int nevts = run->GetNumberOfEvent();
  DMSG(1, "G4VUserMPIrunMerger::Send() : Sending a G4run (" << run << ") with " << nevts
                                                            << " events to: " << destination);
  input_userdata.clear();
  Pack();  // User code
  InputUserData(&nevts, MPI_INT, 1);

  DestroyBuffer();
  G4int newbuffsize = 0;
  G4int size_of_type;

  for (const const_registered_data& el : input_userdata) {
    MPI_Type_size(el.dt, &size_of_type);
    newbuffsize += (size_of_type * el.count);
  }
  char* buffer = new char[newbuffsize];
  // Avoid complains from valgrind (i'm not really sure why this is needed, but, beside the
  // small cpu penalty, we can live with that).)
  std::fill(buffer, buffer + newbuffsize, 0);
  ownsBuffer = true;
  SetupOutputBuffer(buffer, newbuffsize, 0);
  DMSG(3, "Buffer size: " << newbuffsize << " bytes at: " << (void*)outputBuffer);

  // Now userdata contains all data to be send, do the real packing
  for (const const_registered_data& el : input_userdata) {
#ifdef G4MPI_USE_MPI_PACK_NOT_CONST
    MPI_Pack(const_cast<void*>(el.p_data), el.count, el.dt,
#else
    MPI_Pack(el.p_data, el.count, el.dt,
#endif
             outputBuffer, outputBufferSize, &outputBufferPosition, COMM_G4COMMAND_);
  }
  assert(outputBufferSize == outputBufferPosition);
  MPI_Send(outputBuffer, outputBufferSize, MPI_PACKED, destination,
                       G4MPImanager::kTAG_RUN, COMM_G4COMMAND_);
  bytesSent += outputBufferSize;
  DMSG(2, "G4VUserMPIrunMerger::Send() : Done ");
}

void G4VUserMPIrunMerger::Receive(const unsigned int source)
{
  const MPI_Comm* parentComm = G4MPImanager::GetManager()->GetComm();
  int rank;
  MPI_Comm_rank(*parentComm, &rank);
  DMSG(1, "G4VUserMPIrunMerger::Receive(...) , this rank : " << rank << " and receiving from : " << source);
  // DestroyBuffer();
  // Receive from all but one
  // for (G4int rank = 0; rank < commSize-1; ++rank)
  //{
  MPI_Status status;
  MPI_Probe(source, G4MPImanager::kTAG_RUN, COMM_G4COMMAND_, &status);
  // const G4int source = status.Get_source();
  int newbuffsize1;
  MPI_Get_count(&status, MPI_PACKED, &newbuffsize1);
  const G4int newbuffsize = newbuffsize1;
  DMSG(2, "Preparing to receive buffer of size: " << newbuffsize);
  char* buffer = outputBuffer;
  if (newbuffsize > outputBufferSize) {
    DMSG(3, "New larger buffer expected, resize");
    // New larger buffer incoming, recreate buffer
    delete[] outputBuffer;
    buffer = new char[newbuffsize];
    // Avoid complains from valgrind (i'm not really sure why this is needed, but, beside the
    // small cpu penalty, we can live with that).)
    std::fill(buffer, buffer + newbuffsize, 0);
    ownsBuffer = true;
  }
  SetupOutputBuffer(buffer, newbuffsize, 0);
  MPI_Recv(buffer, newbuffsize, MPI_PACKED, source, G4MPImanager::kTAG_RUN, COMM_G4COMMAND_, &status);
  DMSG(3, "Buffer Size: " << outputBufferSize << " bytes at: " << (void*)outputBuffer);
  output_userdata.clear();
  // User code, if implemented will return the concrete G4Run class
  G4Run* aNewRun = UnPack();
  if (aNewRun == nullptr) aNewRun = new G4Run;
  // Add number of events counter
  G4int nevets = 0;
  OutputUserData(&nevets, MPI_INT, 1);
  // now userdata contains all data references, do the real unpacking
  for (const registered_data& el : output_userdata) {
    MPI_Unpack(outputBuffer, outputBufferSize, &outputBufferPosition, el.p_data, el.count, el.dt,
               COMM_G4COMMAND_);
  }
  for (G4int i = 0; i < nevets; ++i)
    aNewRun->RecordEvent(nullptr);

  // Now merge received MPI run with global one
  DMSG(2, "Before G4Run::Merge : " << run->GetNumberOfEvent());
  run->Merge(aNewRun);
  DMSG(2, "After G4Run::Merge : " << run->GetNumberOfEvent());
  delete aNewRun;
  //}
}

void G4VUserMPIrunMerger::Merge()
{
  // G4cout << "G4VUserMPIrunMerger::Merge called" << G4endl;

  DMSG(0, "G4VUserMPIrunMerger::Merge called");
  const MPI_Comm* parentComm = G4MPImanager::GetManager()->GetComm();
  int rank;
  MPI_Comm_rank(*parentComm, &rank);
  const unsigned int myrank = rank;
  commSize = G4MPImanager::GetManager()->GetActiveSize();
  // do not include extra worker in this communication

  if (commSize == 1) {
    DMSG(1, "Comm world size is 1, nothing to do");
    return;
  }
  MPI_Comm_dup(*parentComm, &COMM_G4COMMAND_);
  bytesSent = 0;
  const G4double sttime = MPI_Wtime();

  // Use G4MPIutils to optimize communications between ranks
  typedef std::function<void(unsigned int)> handler_t;
  using std::placeholders::_1;
  handler_t sender = std::bind(&G4VUserMPIrunMerger::Send, this, _1);
  handler_t receiver = std::bind(&G4VUserMPIrunMerger::Receive, this, _1);
  std::function<void(void)> barrier = std::bind(MPI_Barrier, COMM_G4COMMAND_);
  // G4cout << "go to  G4mpi::Merge" << G4endl;
  G4mpi::Merge(sender, receiver, barrier, commSize, myrank);

  // OLD Style p2p communications
  /*
    if ( myrank != destinationRank ) {
        DMSG(0,"Comm world size: "<<commSize<<" this rank is: "
             <<myrank<<" sending to rank "<<destinationRank);
        Send(destinationRank);
    } else {
        DMSG(1,"Comm world size: "<<commSize<<" this rank is: "
            <<myrank<<" receiving. ");
        for ( unsigned int i = 0 ; i<commSize ; ++i) {
            if ( i != myrank ) Receive(i);
        }
    }
  */
  const G4double elapsed = MPI_Wtime() - sttime;
  long total = 0;
  MPI_Reduce(&bytesSent, &total, 1, MPI_LONG, MPI_SUM, destinationRank, COMM_G4COMMAND_);
  if (verbose > 0 && myrank == destinationRank) {
    // Collect from ranks how much data was sent around
    G4cout << "G4VUserMPIrunMerger::Merge() - data transfer performances: "
           << double(total) / 1000. / elapsed << " kB/s"
           << " (Total Data Transfer= " << double(total) / 1000 << " kB in " << elapsed << " s)."
           << G4endl;
  }

  MPI_Comm_free(&COMM_G4COMMAND_);
  DMSG(0, "G4VUserMPIrunMerger::Merge done");
}
