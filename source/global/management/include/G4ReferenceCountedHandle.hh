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
// G4ReferenceCountedHandle
//
// Class description:
//
// A class to provide reference counting mechanism.
// It is a templated class, acting as a smart pointer,
// wrapping the type to be counted. It performs the reference counting
// during the life-time of the counted object. When its count reaches zero
// the counted object is destroyed by explicit call to its destructor.
// This class provides overloaded operators *() and ->() to allow similar
// syntax as for the normal "dumb" pointers.
// The basic rule for the use of this class is that a handle must always
// be exchanged by reference never dinamically allocated (i.e. never
// instantiated using 'new').
// The validity of a smart pointer object can be verified by using the
// operator !() or operator bool(). I.e.:
//    if( !smartPtrObj ) { ... } // Problem! We must initialize it first!
//    else               { ... } // OK!
// Trying to 'delete' a smart pointer object will generate a compilation
// error (since we're dealing with objects, not pointers!).

// Author: Radovan Chytracek, CERN - November 2001
// --------------------------------------------------------------------
#ifndef G4REFERENCECOUNTEDHANDLE_HH
#define G4REFERENCECOUNTEDHANDLE_HH 1

#include "G4Allocator.hh"
#include "G4Types.hh"

template <class X>
class G4CountedObject;

template <class X>
class G4ReferenceCountedHandle
{
 public:
  inline G4ReferenceCountedHandle(X* rep = nullptr);
  // Constructor.

  inline G4ReferenceCountedHandle(const G4ReferenceCountedHandle<X>& right);
  // Copy constructor.

  inline ~G4ReferenceCountedHandle();
  // Destructor.

  inline G4ReferenceCountedHandle<X>& operator=(
    const G4ReferenceCountedHandle<X>& right);
  // Assignment operator by reference.

  inline G4ReferenceCountedHandle<X>& operator=(X* objPtr);
  // Assignment operator by pointer.

  inline unsigned int Count() const;
  // Forward to Counter class.

  inline X* operator->() const;
  // Operator -> allowing the access to counted object.
  // The check for 0-ness is left out for performance reasons,
  // see operator () below.
  // May be called on initialised smart-pointer only!

  inline G4bool operator!() const;
  // Validity test operator.

  inline operator bool() const;
  // Boolean operator.

  inline X* operator()() const;
  // Functor operator (for convenience).

  // There is no provision that this class is subclassed.
  // If it is subclassed & new data members are added then the
  // following "new" & "delete" will fail and give errors.
  //
  inline void* operator new(std::size_t);
  // Operator new defined for G4Allocator.

  inline void operator delete(void* pObj);
  // Operator delete defined for G4Allocator.

 private:
  G4CountedObject<X>* fObj = nullptr;
  // The object subject to reference counting.
};

extern G4GLOB_DLL G4Allocator<G4ReferenceCountedHandle<void>>*& aRCHAllocator();

template <class X>
class G4CountedObject
{
  friend class G4ReferenceCountedHandle<X>;

 public:
  G4CountedObject(X* pObj = nullptr);
  // Constructor.

  ~G4CountedObject();
  // Destructor.

  inline void AddRef();
  // Increase the count.

  inline void Release();
  // Decrease the count and if zero destroy itself.

  // There is no provision that this class is subclassed.
  // If it is subclassed & new data members are added then the
  // following "new" & "delete" will fail and give errors.
  //
  inline void* operator new(std::size_t);
  // Operator new defined for G4Allocator.

  inline void operator delete(void* pObj);
  // operator delete defined for G4Allocator.

 private:
  unsigned int fCount = 0;
  // Reference counter.
  X* fRep = nullptr;
  // The counted object.
};

extern G4GLOB_DLL G4Allocator<G4CountedObject<void>>*&
aCountedObjectAllocator();

// --------- G4CountedObject<X> Inline function definitions ---------

template <class X>
G4CountedObject<X>::G4CountedObject(X* pObj)
  : fRep(pObj)
{
  if(pObj != nullptr)
    fCount = 1;
}

template <class X>
G4CountedObject<X>::~G4CountedObject()
{
  delete fRep;
}

template <class X>
void G4CountedObject<X>::AddRef()
{
  ++fCount;
}

template <class X>
void G4CountedObject<X>::Release()
{
  if(--fCount == 0)
    delete this;
}

template <class X>
void* G4CountedObject<X>::operator new(std::size_t)
{
  if(aCountedObjectAllocator() == nullptr)
    aCountedObjectAllocator() = new G4Allocator<G4CountedObject<void>>;
  return ((void*) aCountedObjectAllocator()->MallocSingle());
}

template <class X>
void G4CountedObject<X>::operator delete(void* pObj)
{
  aCountedObjectAllocator()->FreeSingle((G4CountedObject<void>*) pObj);
}

// --------- G4ReferenceCountedHandle<X> Inline function definitions ---------

template <class X>
G4ReferenceCountedHandle<X>::G4ReferenceCountedHandle(X* rep)
{
  if(rep != nullptr)
    fObj = new G4CountedObject<X>(rep);
}

template <class X>
G4ReferenceCountedHandle<X>::G4ReferenceCountedHandle(
  const G4ReferenceCountedHandle<X>& right)
  : fObj(right.fObj)
{
  fObj->AddRef();
}

template <class X>
G4ReferenceCountedHandle<X>::~G4ReferenceCountedHandle()
{
  if(fObj != nullptr)
    fObj->Release();
}

template <class X>
G4ReferenceCountedHandle<X>& G4ReferenceCountedHandle<X>::operator=(
  const G4ReferenceCountedHandle<X>& right)
{
  if(fObj != right.fObj)
  {
    if(fObj != nullptr)
      fObj->Release();
    this->fObj = right.fObj;
    fObj->AddRef();
  }
  return *this;
}

template <class X>
G4ReferenceCountedHandle<X>& G4ReferenceCountedHandle<X>::operator=(X* objPtr)
{
  if(fObj != nullptr)
    fObj->Release();
  this->fObj = new G4CountedObject<X>(objPtr);
  return *this;
}

template <class X>
unsigned int G4ReferenceCountedHandle<X>::Count() const
{
  return ((fObj != nullptr) ? fObj->fCount : 0);
}

template <class X>
X* G4ReferenceCountedHandle<X>::operator->() const
{
  return ((fObj != nullptr) ? fObj->fRep : 0);
}

template <class X>
G4bool G4ReferenceCountedHandle<X>::operator!() const
{
  return fObj == nullptr;
}

template <class X>
G4ReferenceCountedHandle<X>::operator bool() const
{
  return fObj != nullptr;
}

template <class X>
X* G4ReferenceCountedHandle<X>::operator()() const
{
  return ((fObj != nullptr) ? fObj->fRep : nullptr);
}

template <class X>
void* G4ReferenceCountedHandle<X>::operator new(std::size_t)
{
  if(aRCHAllocator() == nullptr)
    aRCHAllocator() = new G4Allocator<G4ReferenceCountedHandle<void>>;
  return ((void*) aRCHAllocator()->MallocSingle());
}

template <class X>
void G4ReferenceCountedHandle<X>::operator delete(void* pObj)
{
  aRCHAllocator()->FreeSingle((G4ReferenceCountedHandle<void>*) pObj);
}

#endif
