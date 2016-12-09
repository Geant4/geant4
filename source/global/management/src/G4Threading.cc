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
// $Id:$
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
// G4Threading.cc
//
// ---------------------------------------------------------------
// Author: Andrea Dotti (15 Feb 2013): First Implementation
// ---------------------------------------------------------------

#include "G4Threading.hh"
#include "G4AutoDelete.hh"
#include "G4AutoLock.hh"
#include "globals.hh"

#if defined (WIN32)
   #include <Windows.h>
#else
   #include <unistd.h>
   #include <sys/types.h>
   #include <sys/syscall.h>
#endif

#if defined(G4MULTITHREADED)

#include <atomic>

namespace
{
   G4ThreadLocal G4int G4ThreadID = G4Threading::MASTER_ID;
   G4bool isMTAppType = false;
}

G4Pid_t G4Threading::G4GetPidId()
{ // In multithreaded mode return Thread ID
   #if defined(__MACH__)
   #if MAC_OS_X_VERSION_MAX_ALLOWED >= MAC_OS_X_VERSION_10_12
     uint64_t tid64;
     pthread_threadid_np(NULL, &tid64);
     return (pid_t)tid64;
   #else
     return syscall(SYS_thread_selfid);
   #endif
   #elif defined(WIN32)
     return GetCurrentThreadId();
   #else
     return syscall(SYS_gettid);
   #endif
}

G4int G4Threading::G4GetNumberOfCores()
{
   #if defined(WIN32)
     SYSTEM_INFO sysinfo;
     GetSystemInfo( &sysinfo );
     return static_cast<G4int>( sysinfo.dwNumberOfProcessors );
   #else
     return static_cast<G4int>(sysconf( _SC_NPROCESSORS_ONLN ));
   #endif
}

void G4Threading::G4SetThreadId(G4int value ) { G4ThreadID = value; }
G4int G4Threading::G4GetThreadId() { return G4ThreadID; }
G4bool G4Threading::IsWorkerThread() { return (G4ThreadID>=0); }
G4bool G4Threading::IsMasterThread() { return (G4ThreadID==MASTER_ID); }

#if defined(WIN32)  // WIN32 stuff needed for MT
DWORD /*WINAPI*/ G4WaitForSingleObjectInf( __in G4Mutex m )
 { return WaitForSingleObject( m , INFINITE); }
BOOL G4ReleaseMutex( __in G4Mutex m)
 { return ReleaseMutex(m); }
#endif

#if defined(__linux__) || defined(_AIX)
G4bool G4Threading::G4SetPinAffinity(G4int cpu, G4Thread& aT)
{
   cpu_set_t* aset = new cpu_set_t;
   G4AutoDelete::Register(aset);
   CPU_ZERO(aset);
   CPU_SET(cpu,aset);
   return ( pthread_setaffinity_np(aT, sizeof(cpu_set_t), aset) == 0 );
}
#else //Not available for Mac, WIN,...
G4bool G4Threading::G4SetPinAffinity(G4int, G4Thread&)
{
   G4Exception("G4Threading::G4SetPinAffinity()",
          "NotImplemented", JustWarning,
          "Affinity setting not available for this architecture, ignoring...");
   return true;
}
#endif

void G4Threading::SetMultithreadedApplication(G4bool value ) { isMTAppType = value; }
G4bool G4Threading::IsMultithreadedApplication() { return isMTAppType; }

namespace
{
  std::atomic_int numActThreads(0);
}
int G4Threading::WorkerThreadLeavesPool() { return numActThreads--; }
int G4Threading::WorkerThreadJoinsPool() { return numActThreads++;}
G4int G4Threading::GetNumberOfRunningWorkerThreads() { return numActThreads.load(); }

#else  // Sequential mode

G4int fake_mutex_lock_unlock( G4Mutex* ) { return 0; }

G4Pid_t G4Threading::G4GetPidId()
{                    // In sequential mode return Process ID and not Thread ID
    #if defined(WIN32)
    return GetCurrentProcessId();
    #else
    return getpid();
    #endif
}

G4int G4Threading::G4GetNumberOfCores() { return 1; }
G4int G4Threading::G4GetThreadId() { return G4Threading::SEQUENTIAL_ID; }
G4bool G4Threading::IsWorkerThread() { return false; }
G4bool G4Threading::IsMasterThread() { return true; }
void G4Threading::G4SetThreadId(G4int) {}

G4bool G4Threading::G4SetPinAffinity(G4int,G4Thread&) { return true;}

void G4Threading::SetMultithreadedApplication(G4bool) {}
G4bool G4Threading::IsMultithreadedApplication() { return false; }
int G4Threading::WorkerThreadLeavesPool() { return 0; }
int G4Threading::WorkerThreadJoinsPool() { return 0; }
G4int G4Threading::GetNumberOfRunningWorkerThreads() { return 0; }

#endif
