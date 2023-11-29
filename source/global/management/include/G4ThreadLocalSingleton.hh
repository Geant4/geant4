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
// G4ThreadLocalSingleton
//
// Class description:
//
// This class implements a thread-private "singleton". Being thread
// private the singleton is not a singleton in the term, but a different
// instance exists for each thread.
// This class is a wrapper around the real object that we need to
// make singleton.
//
// Limitation:
//   The object that is made thread-private singleton should not
//   contain any thread-local data member. Note that in general,
//   if an object is to be thread-private it is unnecessary to mark
//   any data-member as G4ThreadLocal.
//
// Performance issues:
//   This class uses locks and mutexes.
//
// Example:
//   This is the singleton pattern often found in Geant4 (sequential):
//   class G4Class
//   {
//     private:
//       static G4Class* instance;
//        G4Class() { ... }
//     public:
//       static G4Class* GetInstance()
//       {
//         static G4Class theInstance;
//         if ( instance == nullptr ) instance = &theInstance;
//         return instance;
//       }
//   };
//   This is transformed to the following to implement a thread-local
//   singleton:
//   class G4Class
//   {
//     private:
//       static G4ThreadLocal G4Class* instance;
//       G4Class() { ... }
//     public:
//       static G4Class* GetInstance()
//       {
//         if ( instance == nullptr ) instance = new G4Class;
//           return instance;
//       }
//   };
//   Note that this class also has a memory leak.
//
//   This class can be used as follows:
//   class G4Class
//   {
//     friend class G4ThreadLocalSingleton<G4Class>;
//     private:
//       G4Class() { ... }
//     public:
//       static G4Class* GetInstance()
//       {
//         static G4ThreadLocalSingleton<G4Class> instance;
//         return instance.Instance();
//       }
//   };
//   Each thread has its own instance of G4Class.
//   Deletion of G4Class instances is done at end of program.
//   Note the "friend" statement.

// Author: A.Dotti, 28 October 2013
// --------------------------------------------------------------------
#ifndef G4TLSSINGLETON_HH
#define G4TLSSINGLETON_HH 1

#include "G4AutoLock.hh"
#include "G4Cache.hh"
#include "G4Backtrace.hh"
#include "G4Threading.hh"

#include <list>
#include <vector>
#include <functional>

// Forward declaration. See G4AutoDelete.hh
//
namespace G4AutoDelete
{
  template <class T>
  void Register(T*);
}

template <class T>
class G4ThreadLocalSingleton;

// this explicit specialization holds all the callbacks
// to explicitly invoke the auto-deletion
template <>
class G4ThreadLocalSingleton<void>
{
 private:
  using fvector_t = std::vector<std::function<void()>>;

  template <class T>
  friend class G4ThreadLocalSingleton;

  static fvector_t& GetCallbacks();
  static G4Mutex& GetMutex();

 public:
  static void Clear();

  template <typename FuncT>
  static typename fvector_t::iterator Insert(FuncT&& _func)
  {
    G4AutoLock _lk{ GetMutex() };
    return GetCallbacks().emplace(GetCallbacks().end(),
                                  std::forward<FuncT>(_func));
  }
};

template <class T>
class G4ThreadLocalSingleton : private G4Cache<T*>
{
  friend void G4AutoDelete::Register<T>(T*);

 public:
  G4ThreadLocalSingleton();
  // Creates thread-local singleton manager

  ~G4ThreadLocalSingleton() override;

  G4ThreadLocalSingleton(const G4ThreadLocalSingleton&) = delete;
  G4ThreadLocalSingleton(G4ThreadLocalSingleton&&)      = default;

  G4ThreadLocalSingleton& operator=(const G4ThreadLocalSingleton&) = delete;
  G4ThreadLocalSingleton& operator=(G4ThreadLocalSingleton&&) = default;

  T* Instance() const;
  // Returns a pointer to a thread-private instance of T

 private:
  void Register(T* i) const;

  void Clear();

  mutable std::list<T*> instances;
  mutable G4Mutex listm;
};

//=============================================================
// Inline methods implementation
//=============================================================

template <class T>
G4ThreadLocalSingleton<T>::G4ThreadLocalSingleton()
  : G4Cache<T*>()
{
  G4MUTEXINIT(listm);
  G4Cache<T*>::Put(nullptr);
  // Uncomment below to find the origin of where instantiation happened
  /*
  auto bt = G4Backtrace::GetDemangled<4, 1>(
    [](const char* cstr) { return std::string{ cstr }; });
  std::cout << "Backtrace to G4ThreadLocalSingleton<"
            << G4Demangle<T>().c_str() << ">:\n";
  for(auto& itr : bt)
  {
    if(!itr.empty())
      std::cout << "\t" << itr << "\n";
  }
  */
  G4ThreadLocalSingleton<void>::Insert([&]() {
    // printf("Deleting G4ThreadLocalSingletons for type %s ...\n",
    //       G4Demangle<T>().c_str());
    this->Clear();
  });
}

template <class T>
G4ThreadLocalSingleton<T>::~G4ThreadLocalSingleton()
{
  Clear();
  G4MUTEXDESTROY(listm);
}

template <class T>
T* G4ThreadLocalSingleton<T>::Instance() const
{
  T* instance = G4Cache<T*>::Get();
  if(instance == static_cast<T*>(0))
  {
    instance = new T;
    G4Cache<T*>::Put(instance);
    Register(instance);
  }
  return instance;
}

template <class T>
void G4ThreadLocalSingleton<T>::Register(T* i) const
{
  G4AutoLock l(&listm);
  instances.push_back(i);
}

template <class T>
void G4ThreadLocalSingleton<T>::Clear()
{
  if(instances.empty())
    return;
  G4AutoLock l(&listm);
  while(!instances.empty())
  {
    T* thisinst = instances.front();
    instances.pop_front();
    delete thisinst;
  }
}

#endif
