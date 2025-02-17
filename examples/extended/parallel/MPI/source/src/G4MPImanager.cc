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

#include "mpi.h"

#include "G4MPIbatch.hh"
#include "G4MPIextraWorker.hh"
#include "G4MPImessenger.hh"
#include "G4MPIrandomSeedGenerator.hh"
#include "G4MPIsession.hh"
#include "G4MPIstatus.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4UIcommand.hh"
#include "G4UImanager.hh"

#include <getopt.h>
#include <stdio.h>
#include <time.h>
#include <vector>

G4MPImanager* G4MPImanager::g4mpi_ = NULL;

// --------------------------------------------------------------------------
namespace
{

// wrappers for thread functions
void thread_ExecuteThreadCommand(const G4String* command)
{
  G4MPImanager::GetManager()->ExecuteThreadCommand(*command);
}

// --------------------------------------------------------------------------
void Wait(G4int ausec)
{
  struct timespec treq, trem;
  treq.tv_sec = 0;
  treq.tv_nsec = ausec * 1000;

  nanosleep(&treq, &trem);
}

}  // end of namespace

// --------------------------------------------------------------------------
G4MPImanager::G4MPImanager(int nof_extra_workers)
  : verbose_(0),
    COMM_G4COMMAND_(MPI_COMM_NULL),
    processing_comm_(MPI_COMM_NULL),
    collecting_comm_(MPI_COMM_NULL),
    all_comm_(MPI_COMM_NULL),
    qfcout_(false),
    qinitmacro_(false),
    qbatchmode_(false),
    thread_id_(0),
    master_weight_(1.),
    nof_extra_workers_(nof_extra_workers)
{
  // MPI::Init();
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &provided_);
  if (provided_ < MPI_THREAD_SERIALIZED) {
    std::cerr << "Warning: MPI does not provide the requested threading support!" << std::endl;
  }
  Initialize();
}

// --------------------------------------------------------------------------
G4MPImanager::G4MPImanager(int argc, char** argv, int nof_extra_workers)
  : verbose_(0),
    COMM_G4COMMAND_(MPI_COMM_NULL),
    processing_comm_(MPI_COMM_NULL),
    collecting_comm_(MPI_COMM_NULL),
    all_comm_(MPI_COMM_NULL),
    qfcout_(false),
    qinitmacro_(false),
    qbatchmode_(false),
    thread_id_(0),
    master_weight_(1.),
    nof_extra_workers_(nof_extra_workers)
{
  int provided;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED, &provided_);
  if (provided_ < MPI_THREAD_SERIALIZED) {
    std::cerr << "Warning: MPI does not provide the requested threading support!" << std::endl;
  }

  Initialize();
  ParseArguments(argc, argv);
}

// --------------------------------------------------------------------------
G4MPImanager::~G4MPImanager()
{
  if (is_slave_ && qfcout_) fscout_.close();

  delete status_;
  delete messenger_;
  delete session_;

  if (nof_extra_workers_) {
    MPI_Group_free(&world_group_);
    MPI_Group_free(&processing_group_);
    MPI_Group_free(&collecting_group_);
    MPI_Group_free(&all_group_);
    if (processing_comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&processing_comm_);
    }
    if (collecting_comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&collecting_comm_);
    }
    if (all_comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&all_comm_);
    }
  }
  else {
    MPI_Comm_free(&COMM_G4COMMAND_);
  }

  MPI_Finalize();
}

// --------------------------------------------------------------------------
G4MPImanager* G4MPImanager::GetManager()
{
  if (g4mpi_ == NULL) {
    G4Exception("G4MPImanager::GetManager()", "MPI001", FatalException,
                "G4MPImanager is not instantiated.");
  }
  return g4mpi_;
}

// --------------------------------------------------------------------------
void G4MPImanager::SetExtraWorker(G4VMPIextraWorker* extraWorker)
{
  if (!nof_extra_workers_) {
    G4Exception("G4MPImanager::SetExtraWorker()", "MPI001", FatalException,
                "Number of extra workers >0 must be set first.");
  }

  extra_worker_ = extraWorker;
}

// --------------------------------------------------------------------------
void G4MPImanager::Initialize()
{
  // G4cout << "G4MPImanager::Initialize" << G4endl;

  if (g4mpi_ != NULL) {
    G4Exception("G4MPImanager::Initialize()", "MPI002", FatalException,
                "G4MPImanager is already instantiated.");
  }

  g4mpi_ = this;

  // get rank information
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);//world_size_ = MPI::COMM_WORLD.Get_size();
  if (world_size_ - nof_extra_workers_ <= 0) {
    G4Exception("G4MPImanager::SetExtraWorker()", "MPI001", JustWarning,
                "Cannot reserve extra ranks: the MPI size is not sufficient.");
    nof_extra_workers_ = 0;
  }
  size_ = world_size_ - nof_extra_workers_;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_); //rank_ = MPI::COMM_WORLD.Get_rank();
  is_master_ = (rank_ == kRANK_MASTER);
  is_slave_ = (rank_ != kRANK_MASTER);
  is_extra_worker_ = false;

  if (nof_extra_workers_) {
    // G4cout << "Extra workers requested" << G4endl;

    // Define three groups of workers: processing, collecting and all;
    // if no extra workers are declared, all world ranks are processing ranks

    // MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group_);

    // Group 1 - processing ranks
    int* ranks1 = new int[size_];
    for (int i = 0; i < size_; i++)
      ranks1[i] = i;
    // Construct a group containing all of the processing ranks in world_group
    MPI_Group_incl(world_group_, size_, ranks1, &processing_group_);

    // Group 2 - collecting ranks
    int* ranks2 = new int[nof_extra_workers_];
    for (int i = 0; i < nof_extra_workers_; i++)
      ranks2[i] = (world_size_ - nof_extra_workers_) + i;
    // Construct a group containing all of the collecting ranks in world_group
    MPI_Group_incl(world_group_, nof_extra_workers_, ranks2, &collecting_group_);

    // Group 3 - all ranks
    int* ranks3 = new int[world_size_];
    for (int i = 0; i < world_size_; i++)
      ranks3[i] = i;
    // Construct a group containing all of the processing ranks in world_group
    MPI_Group_incl(world_group_, world_size_, ranks3, &all_group_);

    // Create new communicators based on the groups
    MPI_Comm_create_group(MPI_COMM_WORLD, processing_group_, 0, &processing_comm_);
    MPI_Comm_create_group(MPI_COMM_WORLD, collecting_group_, 0, &collecting_comm_);
    MPI_Comm_create_group(MPI_COMM_WORLD, all_group_, 0, &all_comm_);

    // COMM_G4COMMAND_ = processing_comm_ copy
    COMM_G4COMMAND_ = MPI_Comm(processing_comm_);
  }
  else {
    // G4cout << "No extra workers requested" << G4endl;
    // initialize MPI communicator
    MPI_Comm_dup(MPI_COMM_WORLD, &COMM_G4COMMAND_);
  }

  is_extra_worker_ = (collecting_comm_ != MPI_COMM_NULL);

  // new G4MPI stuffs
  messenger_ = new G4MPImessenger();
  messenger_->SetTargetObject(this);
  session_ = new G4MPIsession;
  status_ = new G4MPIstatus;

  if (!is_extra_worker_) {
    // default seed generator is random generator.
    seed_generator_ = new G4MPIrandomSeedGenerator;
    DistributeSeeds();
  }

  // print status of this worker
  // G4cout << this << " world_size_ " << world_size_ << G4endl;
  // G4cout << this << " size_ " << size_ << G4endl;
  // G4cout << this << " nof_extra_workers_  " << nof_extra_workers_  << G4endl;
  // G4cout << this << " is_master_ " << is_master_ << G4endl;
  // G4cout << this << " is_slave_ " << is_slave_ << G4endl;
  // G4cout << this << " is_extra_worker_ " << is_extra_worker_ << G4endl;
  // G4cout << this << " is_processing_worker_ "
  //        <<  (processing_comm_ != MPI_COMM_NULL) << G4endl;
}

// --------------------------------------------------------------------------
void G4MPImanager::ParseArguments(int argc, char** argv)
{
  _options.clear();
  G4int qhelp = 0;
  G4String ofprefix = "mpi";

  G4int optind = -1;

  static struct option long_options[] = {{"help", no_argument, NULL, 'h'}, {"verbose", no_argument, NULL, 'v'}, {"init", required_argument, NULL, 'i'}, {"ofile", optional_argument, NULL, 'o'}, {NULL, 0, NULL, 0}};

  for (int i = 1; i < argc; i++)
  {
    G4String arg = argv[i];

    if (arg == "help")
    {
      qhelp = 1;
    }
    else if (arg == "verbose")
    {
      verbose_ = 1;
    }
    else if (arg == "init")
    {
      qinitmacro_ = true;
      init_file_name_ = optarg;
    }
    else if (arg == "ofile")
    {
      qfcout_ = true;
      if (optarg) ofprefix = optarg;
    }
    else if (arg[0] == 'm')
    {
      G4String sub = arg.substr(0,5);
      if (sub=="macro")
      {
        optind = i;
        qbatchmode_ = true;
      }
    }
    //default:
      //G4cerr << "*** invalid options specified." << G4endl;
      //std::exit(EXIT_FAILURE);
    else
    {
      _options.push_back(arg);
    }
  }

  // show help
  if (qhelp) {
    if (is_master_) ShowHelp();
    MPI_Finalize();
    std::exit(EXIT_SUCCESS);
  }

  // file output
  if (is_slave_ && qfcout_) {
    G4String prefix = ofprefix + ".%03d" + ".cout";
    char str[1024];
    sprintf(str, prefix.c_str(), rank_);
    G4String fname(str);
    fscout_.open(fname.c_str(), std::ios::out);
  }

  // non-option ARGV-elements ...
  if (qbatchmode_) {
    macro_file_name_ = argv[optind];
  }

}

// ====================================================================
void G4MPImanager::UpdateStatus()
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  const G4Run* run = runManager->GetCurrentRun();

  G4int runid, eventid, neventTBP;

  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState g4state = stateManager->GetCurrentState();

  if (run) {
    runid = run->GetRunID();
    neventTBP = run->GetNumberOfEventToBeProcessed();
    eventid = run->GetNumberOfEvent();
    if (g4state == G4State_GeomClosed || g4state == G4State_EventProc) {
      status_->StopTimer();
    }
  }
  else {
    runid = 0;
    eventid = 0;
    neventTBP = 0;
  }

  status_->SetStatus(rank_, runid, neventTBP, eventid, g4state);
}

// --------------------------------------------------------------------------
void G4MPImanager::ShowStatus()
{
  G4int buff[G4MPIstatus::kNSIZE];

  UpdateStatus();
  G4bool gstatus = CheckThreadStatus();

  if (is_master_) {
    status_->Print();  // for maser itself

    G4int nev = status_->GetEventID();
    G4int nevtp = status_->GetNEventToBeProcessed();
    G4double cputime = status_->GetCPUTime();

    // receive from each slave
    for (G4int islave = 1; islave < size_; islave++) {
      MPI_Recv(buff, G4MPIstatus::kNSIZE, MPI_INT, islave, kTAG_G4STATUS, COMM_G4COMMAND_, MPI_STATUS_IGNORE);
      status_->UnPack(buff);
      status_->Print();

      // aggregation
      nev += status_->GetEventID();
      nevtp += status_->GetNEventToBeProcessed();
      cputime += status_->GetCPUTime();
    }

    G4String strStatus;
    if (gstatus) {
      strStatus = "Run";
    }
    else {
      strStatus = "Idle";
    }

    G4cout << "-------------------------------------------------------" << G4endl
           << "* #ranks= " << size_ << "   event= " << nev << "/" << nevtp
           << " state= " << strStatus << " time= " << cputime << "s" << G4endl;
  }
  else {
    status_->Pack(buff);
    MPI_Send(buff, G4MPIstatus::kNSIZE, MPI_INT, kRANK_MASTER, kTAG_G4STATUS, COMM_G4COMMAND_);
  }
}

// ====================================================================
void G4MPImanager::DistributeSeeds()
{
  // Do nothing if not processing worker
  if (is_extra_worker_) return;

  std::vector<G4long> seed_list = seed_generator_->GetSeedList();
  G4Random::setTheSeed(seed_list[rank_]);
}

// --------------------------------------------------------------------------
void G4MPImanager::ShowSeeds()
{
  G4long buff;

  if (is_master_) {
    // print master
    G4cout << "* rank= " << rank_ << " seed= " << G4Random::getTheSeed() << G4endl;
    // receive from each slave
    for (G4int islave = 1; islave < size_; islave++) {
      MPI_Recv(&buff, 1, MPI_LONG, islave, kTAG_G4SEED, COMM_G4COMMAND_, MPI_STATUS_IGNORE);
      G4cout << "* rank= " << islave << " seed= " << buff << G4endl;
    }
  }
  else {  // slaves
    buff = G4Random::getTheSeed();
    MPI_Send(&buff, 1, MPI_LONG, kRANK_MASTER, kTAG_G4SEED, COMM_G4COMMAND_);
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::SetSeed(G4int inode, G4long seed)
{
  if (rank_ == inode) {
    CLHEP::HepRandom::setTheSeed(seed);
  }
}

// ====================================================================
G4bool G4MPImanager::CheckThreadStatus()
{
  unsigned buff;
  unsigned qstatus = 0;

  if (is_master_) {
    qstatus = (thread_id_ != 0);
    // get slave status
    for (G4int islave = 1; islave < size_; islave++) {
      MPI_Request request;
      MPI_Irecv(&buff, 1, MPI_UNSIGNED, islave, kTAG_G4STATUS, COMM_G4COMMAND_, &request);

      int flag = 0;  // Variable to hold the status of the test
      while (!flag) {
        MPI_Test(&request, &flag, MPI_STATUS_IGNORE);  // Check if the operation is complete
        ::Wait(1000);  // Optionally wait, replace with your own wait function if needed
      }

      qstatus |= buff;
    }
  }
  else {
    buff = (thread_id_ != 0);
    MPI_Send(&buff, 1, MPI_UNSIGNED, kRANK_MASTER, kTAG_G4STATUS, COMM_G4COMMAND_);

  }

  // broadcast
  buff = qstatus;  // for master
  MPI_Bcast(&buff, 1, MPI_UNSIGNED, kRANK_MASTER, COMM_G4COMMAND_);
  qstatus = buff;  // for slave

  if (qstatus != 0)
    return true;
  else
    return false;
}

// --------------------------------------------------------------------------
void G4MPImanager::ExecuteThreadCommand(const G4String& command)
{
  // this method is a thread function.
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int rc = UI->ApplyCommand(command);

  G4int commandStatus = rc - (rc % 100);

  switch (commandStatus) {
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
  if (thread_id_) {
    pthread_join(thread_id_, 0);
    thread_id_ = 0;
  }

  return;
}

// --------------------------------------------------------------------------
void G4MPImanager::ExecuteBeamOnThread(const G4String& command)
{
  G4bool threadStatus = CheckThreadStatus();

  if (threadStatus) {
    if (is_master_) {
      G4cout << "G4MPIsession:: beamOn is still running." << G4endl;
    }
  }
  else {  // ok
    static G4String cmdstr;
    cmdstr = command;
    G4int rc = pthread_create(&thread_id_, 0, (Func_t)thread_ExecuteThreadCommand, (void*)&cmdstr);
    if (rc != 0)
      G4Exception("G4MPImanager::ExecuteBeamOnThread()", "MPI003", FatalException,
                  "Failed to create a beamOn thread.");
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::JoinBeamOnThread()
{
  if (thread_id_) {
    pthread_join(thread_id_, 0);
    thread_id_ = 0;
  }
}

// ====================================================================
G4String G4MPImanager::BcastCommand(const G4String& command)
{
  // Do nothing if not processing worker
  if (is_extra_worker_) return G4String("exit");

  enum
  {
    kBUFF_SIZE = 512
  };
  static char sbuff[kBUFF_SIZE];
  command.copy(sbuff, kBUFF_SIZE);
  G4int len = command.size();
  sbuff[len] = '\0';  // no boundary check

  // "command" is not yet fixed in slaves at this time.

  // waiting message exhausts CPU in LAM!
  // COMM_G4COMMAND_->Bcast(sbuff, ssize, MPI::CHAR, RANK_MASTER);

  // another implementation
  if (is_master_) {
    for (G4int islave = 1; islave < size_; islave++) {
      MPI_Send(sbuff, kBUFF_SIZE, MPI_CHAR, islave, kTAG_G4COMMAND, COMM_G4COMMAND_);
    }
  }
  else {
    // try non-blocking receive
    MPI_Request request;
    MPI_Irecv(sbuff, kBUFF_SIZE, MPI_CHAR, kRANK_MASTER, kTAG_G4COMMAND, COMM_G4COMMAND_, &request);

    int flag = 0;  // Variable to hold the status of the test
    while (!flag) {
        MPI_Test(&request, &flag, MPI_STATUS_IGNORE);  // Check if the non-blocking receive is done
        ::Wait(1000);  // Optionally wait for a bit before checking again
    }
  }

  return G4String(sbuff);
}

// ====================================================================
void G4MPImanager::ExecuteMacroFile(const G4String& fname, G4bool qbatch)
{
  G4bool currentmode = qbatchmode_;
  qbatchmode_ = true;
  G4MPIbatch* batchSession = new G4MPIbatch(fname, qbatch);
  batchSession->SessionStart();
  delete batchSession;
  qbatchmode_ = currentmode;
}

// --------------------------------------------------------------------------
void G4MPImanager::BeamOn(G4int nevent, G4bool qdivide)
{
  // G4cout << "G4MPImanager::BeamOn " << nevent << G4endl;

#ifndef G4MULTITHREADED
  G4RunManager* runManager = G4RunManager::GetRunManager();
#endif

  if (qdivide) {  // events are divided
    G4double ntot = master_weight_ + size_ - 1.;
    G4int nproc = G4int(nevent / ntot);
    G4int nproc0 = nevent - nproc * (size_ - 1);
    fevents_in_master = nproc0;
    fevents_in_slave = nproc;

    if (verbose_ > 0 && is_master_) {
      G4cout << "#events in master=" << nproc0 << " / "
             << "#events in slave=" << nproc << G4endl;
    }

    status_->StartTimer();  // start timer

#ifdef G4MULTITHREADED
    G4String str_nevt;
    if (is_master_)
      str_nevt = G4UIcommand::ConvertToString(nproc0);
    else
      str_nevt = G4UIcommand::ConvertToString(nproc);
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/beamOn " + str_nevt);
#else
    if (is_master_)
      runManager->BeamOn(nproc0);
    else
      runManager->BeamOn(nproc);
#endif

    status_->StopTimer();  // stop timer
  }
  else {  // same events are generated in each node (for test use)
    if (verbose_ > 0 && is_master_) {
      G4cout << "#events in master=" << nevent << " / "
             << "#events in slave=" << nevent << G4endl;
    }
    status_->StartTimer();  // start timer

#ifdef G4MULTITHREADED
    G4String str_nevt = G4UIcommand::ConvertToString(nevent);
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand("/run/beamOn " + str_nevt);
#else
    runManager->BeamOn(nevent);
#endif

    status_->StopTimer();  // stop timer
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::WaitBeamOn()
{
  // G4cout << "G4MPImanager::WaitBeamOn" << G4endl;

  // Extra worker
  if (is_extra_worker_) {
    if (extra_worker_) {
      G4cout << "Calling extra_worker " << G4endl;
      extra_worker_->BeamOn();
    }
    else {
      G4cout << " !!!! extra_worker_ is not defined " << G4endl;
    }
    return;
  }

  G4int buff = 0;
  if (qbatchmode_) {  // valid only in batch mode
    if (is_master_) {
      // receive from each slave
      for (G4int islave = 1; islave < size_; islave++) {
        // G4cout << "calling Irecv for islave " << islave << G4endl;
        MPI_Request request;
        MPI_Irecv(&buff, 1, MPI_INT, islave, kTAG_G4STATUS, COMM_G4COMMAND_, &request);

        int flag = 0;  // Variable to hold the status of the test
        while (!flag) {
            MPI_Test(&request, &flag, MPI_STATUS_IGNORE);  // Check if the non-blocking receive is done
            ::Wait(1000);  // Optionally wait for a bit before checking again
        }
      }
    }
    else {
      buff = 1;
      // G4cout << "calling send for i " << kRANK_MASTER << G4endl;
      MPI_Send(&buff, 1, MPI_INT, kRANK_MASTER, kTAG_G4STATUS, COMM_G4COMMAND_);
    }
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::Print(const G4String& message)
{
  if (is_master_) {
    std::cout << message << std::flush;
  }
  else {
    if (qfcout_) {  // output to a file
      fscout_ << message << std::flush;
    }
    else {  // output to stdout
      std::cout << rank_ << ":" << message << std::flush;
    }
  }
}

// --------------------------------------------------------------------------
void G4MPImanager::ShowHelp() const
{
  if (is_slave_) return;

  G4cout << "Geant4 MPI interface" << G4endl;
  G4cout << "usage:" << G4endl;
  G4cout << "<app> [options] [macro file]" << G4endl << G4endl;
  G4cout << "   -h, --help              show this message." << G4endl;
  G4cout << "   -v, --verbose           show verbose message" << G4endl;
  G4cout << "   -i, --init=FNAME        set an init macro file" << G4endl;
  G4cout << "   -o, --ofile[=FNAME]     set slave output to a flie" << G4endl;
  G4cout << G4endl;
}
