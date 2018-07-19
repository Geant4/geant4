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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLAllocationPool.hh
 * \brief Singleton for recycling allocation of instances of a given class
 *
 * \date 2nd October 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLALLOCATIONPOOL_HH
#define G4INCLALLOCATIONPOOL_HH

#if defined(INCL_USE_ALLOCATION_POOL) || defined(INCLXX_IN_GEANT4_MODE)

#include <stack>
#include <new>
#include <cstddef>

namespace G4INCL {

  template<typename T>
    class AllocationPool {
      public:
        static AllocationPool &getInstance() {
          if(!theInstance)
            theInstance = new AllocationPool<T>;
          return *theInstance;
        }

        T *getObject() {
          if(theStack.empty())
            return static_cast<T*>(::operator new(sizeof(T)));
          else {
            T *t = theStack.top();
            theStack.pop();
            return t;
          }
        }

        void recycleObject(T *t) {
          theStack.push(t);
        }

        void clear() {
          while(!theStack.empty()) { /* Loop checking, 10.07.2015, D.Mancusi */
            ::operator delete(theStack.top());
            theStack.pop();
          }
        }

      protected:
        AllocationPool() {}
        virtual ~AllocationPool() {
          clear();
        }

        static G4ThreadLocal AllocationPool *theInstance;

        std::stack<T*> theStack;

    };

  template<typename T>
    G4ThreadLocal AllocationPool<T> *AllocationPool<T>::theInstance = 0;

}

#define INCL_DECLARE_ALLOCATION_POOL(T) \
  public: \
    static void *operator new(size_t /* s */) { \
      ::G4INCL::AllocationPool<T> &allocator = ::G4INCL::AllocationPool<T>::getInstance(); \
      return allocator.getObject(); \
    } \
    static void operator delete(void *a, size_t /* s */) { \
      ::G4INCL::AllocationPool<T> &allocator = ::G4INCL::AllocationPool<T>::getInstance(); \
      allocator.recycleObject(static_cast<T *>(a)); \
    }

#else // defined(INCL_USE_ALLOCATION_POOL) || defined(INCLXX_IN_GEANT4_MODE)
#define INCL_DECLARE_ALLOCATION_POOL(T)
#endif

#endif // G4INCLALLOCATIONPOOL_HH
