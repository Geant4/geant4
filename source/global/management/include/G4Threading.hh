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
// $Id$
//
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// Class Description:
//
// This file defines types and macros used to expose Geant4 threading model.

// ---------------------------------------------------------------
// Author: Andrea Dotti (15 Feb 2013): First Implementation
// ---------------------------------------------------------------
#ifndef G4Threading_hh
#define G4Threading_hh

#include "G4Types.hh"

// Macro to put current thread to sleep
//
#if defined(WIN32)
#define G4THREADSLEEP( tick ) { Sleep(tick); }
#else
#include <unistd.h>    // needed for sleep()
#define G4THREADSLEEP( tick ) { sleep(tick); }
#endif

#if defined(G4MULTITHREADED)
  //===============================
  // Multi-threaded build
  //===============================
#if ( defined(__MACH__) && defined(__clang__) && defined(__x86_64__) ) || \
    ( defined(__MACH__) && defined(__GNUC__) && (__GNUC__>=4 && __GNUC_MINOR__>=7 || __GNUC__>=5) ) || \
    defined(__linux__) || defined(_AIX)
    //
    // Multi-threaded build: for POSIX systems
    //
    #include <pthread.h>
    #if defined(__MACH__)  // needed only for MacOSX for definition of pid_t
      #include <sys/types.h>
    #endif

    typedef pthread_mutex_t G4Mutex;
    typedef pthread_t G4Thread;

    // G4Mutex initializer macro
    //
    #define G4MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER

    // Lock/unlock a G4Mutex function name
    //
    #define G4MUTEXLOCK pthread_mutex_lock
    #define G4MUTEXUNLOCK pthread_mutex_unlock

    // Macro to initialize a Mutex
    //
    #define G4MUTEXINIT(mutex) pthread_mutex_init( &mutex , NULL);
    #define G4MUTEXDESTROY(mutex) pthread_mutex_destroy( &mutex );

    // Macro to create a G4Thread object
    //
    #define G4THREADCREATE( worker , func , arg )  { \
        pthread_attr_t attr; \
        pthread_attr_init(&attr); \
        pthread_attr_setstacksize(&attr,16*1024*1024); \
        pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE); \
        pthread_create( worker, &attr, func , arg ); \
    }

    // Macro to join thread
    //
    #define G4THREADJOIN( worker ) pthread_join( worker , NULL)

    // Macro to retrieve caller thread
    //
    #define G4THREADSELF pthread_self

    // Some useful types
    //
    typedef void* G4ThreadFunReturnType;
    typedef void* G4ThreadFunArgType;
    typedef G4int (*thread_lock)(G4Mutex*);
    typedef G4int (*thread_unlock)(G4Mutex*);

    typedef pid_t G4Pid_t;

    // Conditions
    //
    // See G4MTRunManager for example on how to use these
    // This complication is needed to be portable with WIN32
    // Note that WIN32 requires an additional initialization step.
    // See example code
    //
    typedef pthread_cond_t G4Condition;
    #define G4CONDITION_INITIALIZER PTHREAD_COND_INITIALIZER

    #define G4CONDITIONWAIT( cond, mutex ) pthread_cond_wait( cond , mutex );
    #define G4CONDITIONBROADCAST( cond ) pthread_cond_broadcast( cond );

  #elif defined(WIN32)
    //
    // Multi-threaded build: for Windows systems
    //
    #include "windefs.hh"  // Include 'safe...' <windows.h>
	
    typedef HANDLE G4Mutex;
    typedef HANDLE G4Thread;

    #define G4MUTEX_INITIALIZER CreateMutex(NULL,FALSE,NULL)
    DWORD /*WINAPI*/ G4WaitForSingleObjectInf( __in G4Mutex m );
    #define G4MUTEXLOCK G4WaitForSingleObjectInf

    // #define G4MUTEXINIT(mutex) InitializeCriticalSection( &mutex );
    #define G4MUTEXINIT(mutex);
    #define G4MUTEXDESTROY(mutex);

    // Not clear why following two lines are needed...
    //
    BOOL G4ReleaseMutex( __in G4Mutex m);
    #define G4MUTEXUNLOCK G4ReleaseMutex

    #define G4THREADCREATE( worker, func, arg ) { *worker = CreateThread( NULL, 16*1024*1024 , func , arg , 0 , NULL ); }
    #define G4THREADJOIN( worker ) WaitForSingleObject( worker , INFINITE);
    #define G4THREADSELF GetCurrentThreadId
    #define G4ThreadFunReturnType DWORD WINAPI
    typedef LPVOID G4ThreadFunArgType;
    typedef DWORD (*thread_lock)(G4Mutex);
    typedef BOOL (*thread_unlock)(G4Mutex);
    typedef DWORD G4Pid_t;

    // Conditions
    //
    typedef CONDITION_VARIABLE G4Condition;
    #define G4CONDITION_INITIALIZER CONDITION_VARIABLE_INIT

    #define G4CONDITIONWAIT( cond , criticalsectionmutex ) SleepConditionVariableCS( cond, criticalsectionmutex , INFINITE );
    #define G4CONDITIONBROADCAST( cond ) WakeAllConditionVariable( cond );

  #else

    #error "No Threading model technology supported for this platform. Use sequential build !"

  #endif

#else  
  //==========================================
  // G4MULTITHREADED is OFF - Sequential build
  //==========================================
  typedef G4int G4Mutex;
  typedef G4int G4Thread;
  #define G4MUTEX_INITIALIZER 1
  G4int fake_mutex_lock_unlock( G4Mutex* );// { return 0; }
  #define G4MUTEXINIT(mutex) ;;
  #define G4MUTEXDESTROY(mutex) ;;
  #define G4MUTEXLOCK fake_mutex_lock_unlock
  #define G4MUTEXUNLOCK fake_mutex_lock_unlock
  #define G4THREADCREATE( worker , func , arg ) ;;
  #define G4THREADJOIN( worker ) ;;
  #define G4THREADSELF( nothing ) G4Thread(nothing); 
  typedef void* G4ThreadFunReturnType;
  typedef void* G4ThreadFunArgType;
  typedef G4int (*thread_lock)(G4Mutex*);
  typedef G4int (*thread_unlock)(G4Mutex*);
  typedef G4int G4Pid_t;
  typedef G4int G4Condition;
  #define G4CONDITION_INITIALIZER 1
  #define G4CONDITIONWAIT( cond, mutex ) { ++(*cond); ++(*mutex); }
  #define G4CONDITIONBROADCAST( cond ) { ++(*cond); }

#endif //G4MULTITHREADING

namespace G4Threading
{
  enum {
        SEQUENTIAL_ID = -2,
        MASTER_ID = -1,
        WORKER_ID = 0,
        GENERICTHREAD_ID = -1000
    };
  G4Pid_t G4GetPidId();
  G4int G4GetNumberOfCores();
  G4int G4GetThreadId();
  G4bool IsWorkerThread();
  G4bool IsMasterThread();
  void G4SetThreadId( G4int aNewValue );
  G4bool G4SetPinAffinity( G4int idx , G4Thread& at);
  void SetMultithreadedApplication(G4bool value);
  G4bool IsMultithreadedApplication();
  int WorkerThreadLeavesPool();
  int WorkerThreadJoinsPool();
  G4int GetNumberOfRunningWorkerThreads();
}

#endif //G4Threading_hh
