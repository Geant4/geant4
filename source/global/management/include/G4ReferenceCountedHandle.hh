//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ReferenceCountedHandle.hh,v 1.13 2001-11-07 00:32:06 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4ReferenceCountedHandle
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

// Author:      Radovan Chytracek, CERN  (Radovan.Chytracek@cern.ch)
// Version:     3.0
// Date:        November 2001
// ----------------------------------------------------------------------
#ifndef _G4REFERENCECOUNTEDHANDLE_H_
#define _G4REFERENCECOUNTEDHANDLE_H_ 1

#include "G4Allocator.hh"

template <class X> class CountedObject;

template <class X>
class G4ReferenceCountedHandle
{

public: // with description

  inline G4ReferenceCountedHandle( X* rep = 0 );
    // Constructor.
  
  inline G4ReferenceCountedHandle( const G4ReferenceCountedHandle<X>& right );
    // Copy constructor.
  
  inline ~G4ReferenceCountedHandle();
    // Destructor.
  
  inline G4ReferenceCountedHandle<X>&
    operator =( const G4ReferenceCountedHandle<X>& right );
    // Assignment operator by reference.
  
  inline G4ReferenceCountedHandle<X>& operator =( X* objPtr );
    // Assignment operator by pointer.
  
  inline unsigned int Count() const;
    // Forward to Counter class.
  
  inline X* operator ->() const;
    // Operator -> allowing the access to counted object.
    // The check for 0-ness is left out for performance reasons,
    // see operator () below.
    // May be called on initialised smart-pointer only!
  
  inline G4bool operator !() const;
    // Validity test operator.
  
  inline operator bool() const;
    // Boolean operator.
  
  inline X* operator ()() const;
    // Functor operator (for convenience).
  
  // There is no provision that this class is subclassed.
  // If it is subclassed & new data members are added then the
  // following "new" & "delete" will fail and give errors. 
  //
  inline void* operator new( size_t );
    // Operator new defined for G4Allocator.
  
  inline void operator delete( void *pObj );
    // Operator delete defined for G4Allocator.
  
  void* operator new( size_t, void *pObj );
    // This is required when this class is used in the context of STL container.
    // On Windows/VC++ and Digital's cxx it will compile with
    // a warning saying something about not existing correspondent delete...

private:

  CountedObject<X>*     fObj;
    // The object subject to reference counting.

  static G4Allocator<G4ReferenceCountedHandle<X> > aRCHAllocator;
};

template <class X>
class CountedObject
{

  friend class G4ReferenceCountedHandle<X>;

public:  // with description

  CountedObject( X* pObj = 0 );
    // Constructor.

  ~CountedObject();
    // Destructor.

  inline void AddRef();
    // Increase the count.

  inline void Release();
    // Decrease the count and if zero destroy itself.
  
  // There is no provision that this class is subclassed.
  // If it is subclassed & new data members are added then the
  // following "new" & "delete" will fail and give errors.
  //
  inline void* operator new( size_t );
    // Operator new defined for G4Allocator.

  inline void operator delete( void *pObj );
    // operator delete defined for G4Allocator.

private:

  unsigned int fCount;
    // Reference counter.
  X* fRep;
    // The counted object.

  static G4Allocator<CountedObject<X> > aCountedObjectAllocator;
};


// --------- CountedObject<X> Inline function definitions ---------

template <class X>
CountedObject<X>::CountedObject( X* pObj )
 : fCount(0), fRep( pObj )
{
    if( pObj != 0 ) {
      fCount = 1;
    }
}

template <class X>
CountedObject<X>::~CountedObject()
{
    delete fRep;
}
    
template <class X>
void CountedObject<X>::AddRef()
{
    ++fCount;
}
    
template <class X>
void CountedObject<X>::Release()
{
    if( --fCount == 0 ) delete this;
}

template <class X>
void* CountedObject<X>::operator new( size_t )
{
    return( (void *)aCountedObjectAllocator.MallocSingle() );
}
    
template <class X>
void CountedObject<X>::operator delete( void *pObj )
{
    aCountedObjectAllocator.FreeSingle( (CountedObject<X>*)pObj );
}


// --------- G4ReferenceCountedHandle<X> Inline function definitions ---------

template <class X>
G4ReferenceCountedHandle<X>::
 G4ReferenceCountedHandle( X* rep )
 : fObj( 0 )
{
  if( rep != 0 ) {
      fObj = new CountedObject<X>( rep );
  }
}

template <class X>
G4ReferenceCountedHandle<X>::
 G4ReferenceCountedHandle( const G4ReferenceCountedHandle<X>& right )
 : fObj( right.fObj )
{
    fObj->AddRef();
}
  
template <class X>
G4ReferenceCountedHandle<X>::~G4ReferenceCountedHandle()
{
    if( fObj ) fObj->Release();
}
  
template <class X>
G4ReferenceCountedHandle<X>& G4ReferenceCountedHandle<X>::
 operator =( const G4ReferenceCountedHandle<X>& right )
{
    if( fObj != right.fObj ) {
      if( fObj )
        fObj->Release();
      this->fObj = right.fObj;
      fObj->AddRef();
    }
    return *this;
}
  
template <class X>
G4ReferenceCountedHandle<X>& G4ReferenceCountedHandle<X>::
 operator =( X* objPtr )
{
    if( fObj )
      fObj->Release();
    this->fObj = new  CountedObject<X>( objPtr );
    return *this;
}
  
template <class X>
unsigned int G4ReferenceCountedHandle<X>::Count() const
{
    return( fObj ? fObj->fCount : 0 );
}
  
template <class X>
X* G4ReferenceCountedHandle<X>::operator ->() const
{
    return( fObj ? fObj->fRep : 0 );
}
  
template <class X>
G4bool G4ReferenceCountedHandle<X>::operator !() const
{
    return( ( !fObj ) ? true : false );
}
  
template <class X>
G4ReferenceCountedHandle<X>::operator bool() const
{
    return( ( fObj ) ? true : false );
}
  
template <class X>
X* G4ReferenceCountedHandle<X>::operator ()() const
{
    return( fObj ? fObj->fRep : 0 );
}
  
template <class X>
void* G4ReferenceCountedHandle<X>::operator new( size_t )
{
    return( (void *)aRCHAllocator.MallocSingle() );
}
  
template <class X>
void G4ReferenceCountedHandle<X>::operator delete( void *pObj )
{
    aRCHAllocator.FreeSingle( (G4ReferenceCountedHandle<X>*)pObj );
}

template <class X>
void* G4ReferenceCountedHandle<X>::operator new( size_t, void *pObj )
{
    return pObj;
}

// ------------------------------------------------------------------

// In order to save the human's typing and brain the macro is provided
// for definition of the allocators for generic type of counted objects.

template <class X>
G4Allocator<CountedObject<X> >
  CountedObject<X>::aCountedObjectAllocator;

template <class X>
G4Allocator<G4ReferenceCountedHandle<X> >
  G4ReferenceCountedHandle<X>::aRCHAllocator;

#endif // _G4REFERENCECOUNTEDHANDLE_H_

