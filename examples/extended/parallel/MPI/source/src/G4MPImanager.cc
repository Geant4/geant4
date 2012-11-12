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
/// @file G4MPImanager.cc
/// @brief MPI manager class

#include "G4MPImanager.hh"
#include "G4MPImessenger.hh"
#include "G4MPIsession.hh"
#include "G4MPIbatch.hh"
#include "G4MPIstatus.hh"
#include "G4MPIrandomSeedGenerator.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4Run.hh"
#include <time.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>

G4MPImanager* G4MPImanager::theManager = 0;

// --------------------------------------------------------------------------
// wrappers for thread functions
static void thread_ExecuteThreadCommand(const G4String* command)
{
  G4MPImanager::GetManager()-> ExecuteThreadCommand(*command);
}

// --------------------------------------------------------------------------
G4MPImanager::G4MPImanager()
  : verbose(0), qfcout(false), qinitmacro(false), qbatchmode(false),
    threadID(0), masterWeight(1.)
{
  //MPI::Init();
  MPI::Init_thread(MPI::THREAD_SERIALIZED);
  Initialize();
}

// --------------------------------------------------------------------------
G4MPImanager::G4MPImanager(int argc, char** argv)
  : verbose(0), qfcout(false), qinitmacro(false), qbatchmode(false),
    threadID(0), masterWeight(1.)
{
  //MPI::Init(argc, argv);
  MPI::Init_thread(argc, argv, MPI::THREAD_SERIALIZED);
  Initialize();
  ParseArguments(argc, argv);
}

// --------------------------------------------------------------------------
G4MPImanager::~G4MPImanager()
{
  if(isSlave && qfcout) fscout.close();

  delete status;
  delete messenger;
  delete session;

  COMM_G4COMMAND.Free();

  MPI::Finalize();

  theManager = 0;
}

// --------------------------------------------------------------------------
G4MPImanager* G4MPImanager::GetManager()
{
  assert( theManager != 0 );
  return theManager;
}

// --------------------------------------------------------------------------
void G4MPImanager::Initialize()
{
  assert( theManager == 0 );

  theManager= this;

  // get rank information
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  isMaster = (rank == RANK_MASTER);
  isSlave = (rank != RANK_MASTER);

  // initialize MPI communicator
  COMM_G4COMMAND = MPI::COMM_WORLD.Dup();

  // new G4MPI stuffs
  messenger = new G4MPImessenger(this);
  session = new G4MPIsession;
  status = new G4MPIstatus;

  // default seed generator is random generator.
  seedGenerator = new G4MPIrandomSeedGenerator;
  DistributeSeeds();
}

// --------------------------------------------------------------------------
void G4MPImanager::ParseArguments(int argc, char** argv)
{
  G4int qhelp = 0;
  G4String ofprefix = "mpi";

  G4int c;
  while (1) {
    G4int option_index= 0;
    static struct option long_options[] = {
      {"help", 0, 0, 0},
      {"verbose", 0, 0, 0},
      {"init", 1, 0, 0},
      {"ofile", 2, 0, 0},
      {0, 0, 0, 0}
    };

    opterr = 0; // suppress message
    c= getopt_long(argc, argv, "hvi:o", long_options, &option_index);
    opterr = 1;

    if(c == -1) break;

    switch (c) {
    case 0:
      switch(option_index) {
      case 0 : // --help
        qhelp = 1;
        break;
      case 1 : // --verbose
        verbose = 1;
        break;
      case 2 : // --init
        qinitmacro = true;
        initFileName = optarg;
        break;
      case 3 : // --ofile
        qfcout = true;
        if(optarg) ofprefix = optarg;
        break;
      }
      break;
    case 'h' :
      qhelp = 1;
      break;
    case 'v' :
      verbose = 1;
      break;
    case 'i' :
      qinitmacro = true;
      initFileName = optarg;
      break;
    case 'o' :
      qfcout = true;
      break;
    default:
      break;
    }
  }

  // show help
  if(qhelp) {
    if(isMaster) ShowHelp();
    MPI::Finalize();
    exit(0);
  }

  // file output
  if(isSlave && qfcout) {
    G4String prefix = ofprefix + ".%03d" + ".cout";
    char str[1024];
    sprintf(str, prefix.c_str(), rank);
    G4String fname(str);
    fscout.open(fname.c_str(), std::ios::out);
  }

  // non-option ARGV-elements ...
  if (optind < argc ) {
    qbatchmode = true;
    macroFileName = argv[optind];
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::Wait(G4int ausec) const
{
  struct timespec treq, trem;
  treq.tv_sec = 0;
  treq.tv_nsec = ausec*1000;

  nanosleep(&treq, &trem);
}

// ====================================================================
void G4MPImanager::UpdateStatus()
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  const G4Run* run = runManager-> GetCurrentRun();

  G4int runid, eventid, neventTBP;

  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState g4state = stateManager-> GetCurrentState();

  if (run) {
    runid = run-> GetRunID();
    neventTBP = run -> GetNumberOfEventToBeProcessed();
    eventid = run-> GetNumberOfEvent();
    if(g4state == G4State_GeomClosed || g4state == G4State_EventProc) {
      status-> StopTimer();
    }
  } else {
    runid = 0;
    eventid = 0;
    neventTBP = 0;
  }

  status-> SetStatus(rank, runid, neventTBP, eventid, g4state);
}

// --------------------------------------------------------------------------
void G4MPImanager::ShowStatus()
{
  G4int buff[G4MPIstatus::NSIZE];

  UpdateStatus();
  G4bool gstatus = CheckThreadStatus();

  if(isMaster) {
    status-> Print();  // for maser itself

    G4int nev = status-> GetEventID();
    G4int nevtp = status-> GetNEventToBeProcessed();
    G4double cputime = status-> GetCPUTime();

    // receive from each slave
    for (G4int islave = 1; islave < size; islave++) {
      COMM_G4COMMAND.Recv(buff, G4MPIstatus::NSIZE, MPI::INT,
                          islave, TAG_G4STATUS);
      status-> UnPack(buff);
      status-> Print();

      // aggregation
      nev += status-> GetEventID();
      nevtp += status-> GetNEventToBeProcessed();
      cputime += status-> GetCPUTime();
    }

    G4String strStatus;
    if(gstatus) {
      strStatus = "Run";
    } else {
      strStatus = "Idle";
    }

    G4cout << "-------------------------------------------------------"
           << G4endl
           << "* #ranks= " << size
           << "   event= " << nev << "/" << nevtp
           << " state= " << strStatus
           << " time= " << cputime << "s"
           << G4endl;
  } else {
    status-> Pack(buff);
    COMM_G4COMMAND.Send(buff, G4MPIstatus::NSIZE, MPI::INT,
                        RANK_MASTER, TAG_G4STATUS);
  }
}

// ====================================================================
void G4MPImanager::DistributeSeeds()
{
  std::vector<G4long> seedList = seedGenerator-> GetSeedList();
  CLHEP::HepRandom::setTheSeed(seedList[rank]);
}

// --------------------------------------------------------------------------
void G4MPImanager::ShowSeeds()
{
  G4long buff;

  if(isMaster) {
    // print master
    G4cout << "* rank= " << rank
           << " seed= " << CLHEP::HepRandom::getTheSeed()
           << G4endl;
    // receive from each slave
    for (G4int islave = 1; islave < size; islave++) {
      COMM_G4COMMAND.Recv(&buff, 1, MPI::LONG, islave, TAG_G4SEED);
      G4cout << "* rank= " << islave
             << " seed= " << buff
             << G4endl;
    }
  } else { // slaves
    buff = CLHEP::HepRandom::getTheSeed();
    COMM_G4COMMAND.Send(&buff, 1, MPI::LONG, RANK_MASTER, TAG_G4SEED);
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::SetSeed(G4int inode, G4long seed)
{
  if(rank == inode) {
    CLHEP::HepRandom::setTheSeed(seed);
  }
}

// ====================================================================
G4bool G4MPImanager::CheckThreadStatus()
{
  unsigned buff;
  G4bool qstatus= false;

  if(isMaster) {
    qstatus = threadID;
    // get slave status
    for (G4int islave = 1; islave < size; islave++) {
      MPI::Request request = COMM_G4COMMAND.Irecv(&buff, 1, MPI::UNSIGNED,
                                                  islave, TAG_G4STATUS);
      MPI::Status status;
      while(! request.Test(status)) {
        Wait(1000);
      }
      qstatus |= buff;
    }
  } else {
    buff = unsigned(threadID);
    COMM_G4COMMAND.Send(&buff, 1, MPI::UNSIGNED, RANK_MASTER, TAG_G4STATUS);
  }

  // broadcast
  buff = qstatus; // for master
  COMM_G4COMMAND.Bcast(&buff, 1, MPI::UNSIGNED, RANK_MASTER);
  qstatus = buff; // for slave

  return qstatus;
}

// --------------------------------------------------------------------------
void G4MPImanager::ExecuteThreadCommand(const G4String& command)
{
  // this method is a thread function.
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int rc = UI-> ApplyCommand(command);

  G4int commandStatus = rc - (rc%100);

  switch(commandStatus) {
  case fCommandSucceeded:
    break;
  case fIllegalApplicationState:
    G4cerr << "illegal application state -- command refused" << G4endl;
    break;
  default:
    G4cerr << "command refused (" << commandStatus << ")" << G4endl;
    break;
  }

  // thread is joined
  if(threadID) {
    pthread_join(threadID, 0);
    threadID = 0;
  }

  return;
}

// --------------------------------------------------------------------------
void G4MPImanager::ExecuteBeamOnThread(const G4String& command)
{
  G4bool threadStatus = CheckThreadStatus();

  if (threadStatus) {
    if(isMaster) {
      G4cout << "G4MPIsession:: beamOn is still running." << G4endl;
    }
  } else { // ok
    static G4String cmdstr;
    cmdstr = command;
    G4int rc = pthread_create(&threadID, 0,
                              (Func_t)thread_ExecuteThreadCommand,
                              (void*)&cmdstr);
    if (rc != 0)
      G4Exception("G4MPImanager::ExecuteBeamOnThread()",
                  "MPI001",
                  FatalException,
                  "Failed to create a beamOn thread.");
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::JoinBeamOnThread()
{
  if(threadID) {
    pthread_join(threadID, 0);
    threadID = 0;
  }
}

// ====================================================================
G4String G4MPImanager::BcastCommand(const G4String& command)
{
  enum { BUFF_SIZE = 512 };
  static char sbuff[BUFF_SIZE];
  command.copy(sbuff,BUFF_SIZE);
  G4int len = command.size();
  sbuff[len] ='\0';   // no boundary check

  // "command" is not yet fixed in slaves at this time.

  // waiting message exhausts CPU in LAM!
  //COMM_G4COMMAND.Bcast(sbuff, ssize, MPI::CHAR, RANK_MASTER);

  // another implementation
  if( isMaster ) {
    for (G4int islave = 1; islave < size; islave++) {
      COMM_G4COMMAND.Send(sbuff, BUFF_SIZE, MPI::CHAR, islave, TAG_G4COMMAND);
    }
  } else {
    // try non-blocking receive
    MPI::Request request= COMM_G4COMMAND.Irecv(sbuff, BUFF_SIZE, MPI::CHAR,
                                               RANK_MASTER, TAG_G4COMMAND);
    // polling...
    MPI::Status status;
    while(! request.Test(status)) {
      Wait(1000);
    }
  }

  return G4String(sbuff);
}

// ====================================================================
void G4MPImanager::ExecuteMacroFile(const G4String& fname, G4bool qbatch)
{
  G4bool currentmode = qbatchmode;
  qbatchmode = true;
  G4MPIbatch* batchSession = new G4MPIbatch(fname, qbatch);
  batchSession-> SessionStart();
  delete batchSession;
  qbatchmode = currentmode;
}

// --------------------------------------------------------------------------
void G4MPImanager::BeamOn(G4int nevent, G4bool qdivide)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();

  if(qdivide) { // events are divided
    G4double ntot = masterWeight+size-1.;
    G4int nproc = G4int(nevent/ntot);
    G4int nproc0 = nevent-nproc*(size-1);

    if(verbose>0 && isMaster) {
      G4cout << "#events in master=" << nproc0 << " / "
             << "#events in slave="  << nproc << G4endl;
    }

    status-> StartTimer(); // start timer
    if(isMaster) runManager-> BeamOn(nproc0);
    else runManager-> BeamOn(nproc);
    status-> StopTimer(); // stop timer

  } else { // same events are generated in each node (for test use)
    if(verbose>0 && isMaster) {
      G4cout << "#events in master=" << nevent << " / "
             << "#events in slave="  << nevent << G4endl;
    }
    status-> StartTimer(); // start timer
    runManager-> BeamOn(nevent);
    status-> StopTimer(); // stop timer
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::WaitBeamOn()
{
  G4int buff = 0;
  if (qbatchmode) {  // valid only in batch mode
    if(isMaster) {
      // receive from each slave
      for (G4int islave = 1; islave < size; islave++) {
        MPI::Request request = COMM_G4COMMAND.Irecv(&buff, 1, MPI::INT,
                                                    islave, TAG_G4STATUS);
        MPI::Status status;
        while(! request.Test(status)) {
           Wait(1000);
        }
      }
    } else {
      buff = 1;
      COMM_G4COMMAND.Send(&buff, 1, MPI::INT, RANK_MASTER, TAG_G4STATUS);
    }
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::Print(const G4String& message)
{
  if(isMaster){
    std::cout << message << std::flush;
  } else {
    if(qfcout) { // output to a file
      fscout << message << std::flush;
    } else { // output to stdout
      std::cout << rank << ":" << message << std::flush;
    }
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::ShowHelp() const
{
  if(isSlave) return;

  G4cout << "Geant4 MPI interface" << G4endl;
  G4cout << "usage:" << G4endl;
  G4cout << "<app> [options] [macro file]"
         << G4endl << G4endl;
  G4cout << "   -h, --help              show this message."
         << G4endl;
  G4cout << "   -v, --verbose           show verbose message"
         << G4endl;
  G4cout << "   -i, --init=FNAME        set an init macro file"
         << G4endl;
  G4cout << "   -o, --ofile[=FNAME]     set slave output to a flie"
         << G4endl;
  G4cout << G4endl;
}
