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

#if defined(G4MULTITHREADED)
  //===============================
  // Multi-threaded build
  //===============================
#if ( defined(__MACH__) && defined(__clang__) && defined(__x86_64__) ) || \
    ( defined(__MACH__) && defined(__GNUC__) && __GNUC__>=4 && __GNUC_MINOR__>=7 ) || \
    defined(__linux__) || defined(_AIX)
    //
    // Multi-threaded build: for POSIX systems
    //
    #include <pthread.h>
    #if defined(__MACH__)
       //needed only for MacOSX for definition of pid_t
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

    // Macro to create a G4Thread object
    //
    #define G4THREADCREATE( worker , func , arg )  { \
        pthread_attr_t attr; \
        pthread_attr_init(&attr); \
        pthread_attr_setstacksize(&attr,16*1024*1024); \
        pthread_create( worker, &attr, func , arg ); \
    }

    // Macro to join a G4Thread
    //
    #define G4THREADJOIN( worker ) pthread_join( worker , NULL)

    // Some useful types
    //
    typedef void* G4ThreadFunReturnType;
    typedef void* G4ThreadFunArgType;
    typedef G4int (*thread_lock)(G4Mutex*);
    typedef G4int (*thread_unlock)(G4Mutex*);

    typedef pid_t G4Pid_t;
  #elif defined(WIN32)
    //
    // Multi-threaded build: for Windows systems
    //
    #include <Windows.h>

    typedef HANDLE G4Mutex;
    typedef HANDLE G4Thread;

    #define G4MUTEX_INITIALIZER CreateMutex(NULL,FALSE,NULL)
    DWORD /*WINAPI*/ G4WaitForSingleObjectInf( __in G4Mutex m );// { return WaitForSingleObject( m , INFINITE); }
    #define G4MUTEXLOCK G4WaitForSingleObjectInf
    //Don't ask me why I need the following two lines...
    BOOL G4ReleaseMutex( __in G4Mutex m);// { return ReleaseMutex(m); }
    #define G4MUTEXUNLOCK G4ReleaseMutex
    #define G4THREADCREATE( worker, func, arg ) { *worker = CreateThread( NULL, 16*1024*1024 , func , arg , 0 , NULL ); }
    #define G4THREADJOIN( worker ) WaitForSingleObject( worker , INFINITE);
    #define G4ThreadFunReturnType DWORD WINAPI
    typedef LPVOID G4ThreadFunArgType;
    typedef DWORD (*thread_lock)(G4Mutex);
    typedef BOOL (*thread_unlock)(G4Mutex);
    typedef DWORD G4Pid_t;
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
  #define G4MUTEXLOCK fake_mutex_lock_unlock
  #define G4MUTEXUNLOCK fake_mutex_lock_unlock
  #define G4THREADCREATE( worker , func , arg ) ;;
  #define G4THREADJOIN( worker ) ;;
  typedef void* G4ThreadFunReturnType;
  typedef void* G4ThreadFunArgType;
  typedef G4int (*thread_lock)(G4Mutex*);
  typedef G4int (*thread_unlock)(G4Mutex*);
  typedef G4int G4Pid_t;
#endif //G4MULTITHREADING

G4Pid_t G4GetPidId();

#endif //G4Threading_hh