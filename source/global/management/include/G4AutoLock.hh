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
// This class provides a mechanism to create a mutex and locks/unlocks it.
// Can be used by applications to implement in a portable way a mutexing logic.
// Usage Example:
//
//      #include "G4Threading.hh"
//      #include "G4AutoLock.hh"
//      /* somehwere */
//      G4Mutex aMutex = G4MUTEX_INITIALIZER;
//      /*
//       somewhere else:
//       The G4AutoLock instance will automatically unlock the mutex when it
//       goes out of scope, lock and unlock method are anyway available for
//       explicit handling of mutex lock. */
//      G4AutoLock l(&aMutex);
//      ProtectedCode();
//      l.unlock(); //explicit unlock
//      UnprotectedCode();
//      l.lock();   //explicit lock
//
// Note that G4AutoLock is defined also for a sequential Geant4 build,
// but has no effect.

// ---------------------------------------------------------------
// Author: Andrea Dotti (15 Feb 2013): First Implementation
// ---------------------------------------------------------------

#ifndef G4AUTOLOCK_HH
#define G4AUTOLOCK_HH

#include "G4Threading.hh"

// Note: Note that G4TemplateAutoLock by itself is not thread-safe and
//       cannot be shared among threads due to the locked switch
//
template<class M, typename L, typename U>
class G4TemplateAutoLock
{
  public:

    G4TemplateAutoLock(M* mtx, L l, U u) : locked(false), _m(mtx), _l(l), _u(u)
    {
        lock();
    }

    virtual ~G4TemplateAutoLock()
    {
        unlock();
    }

    inline void unlock() {
        if ( !locked ) return;
        _u(_m);
        locked = false;
    }

    inline void lock() {
        if ( locked ) return;
        _l(_m);
        locked = true;
    }

  private:

    // Disable copy and assignement operators
    //
    G4TemplateAutoLock( const G4TemplateAutoLock& rhs );
    G4TemplateAutoLock& operator= ( const G4TemplateAutoLock& rhs );

  private:
    G4bool locked;
    M* _m;
    L _l;
    U _u;
};

struct G4ImpMutexAutoLock
  : public G4TemplateAutoLock<G4Mutex,thread_lock,thread_unlock>
{
    G4ImpMutexAutoLock(G4Mutex* mtx)
      : G4TemplateAutoLock<G4Mutex, thread_lock, thread_unlock>
        (mtx, &G4MUTEXLOCK, &G4MUTEXUNLOCK) {}
};
typedef G4ImpMutexAutoLock G4AutoLock;

#endif //G4AUTOLOCK_HH
