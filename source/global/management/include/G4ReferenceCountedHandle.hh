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
// $Id: G4ReferenceCountedHandle.hh,v 1.10 2001-11-03 00:27:28 rado Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4RereferenceCountedHandle
// Class description:
// A class to provide reference counting mechanism for a class which
// is designed without this support. It is a templated class, a smart pointer,
// wrapping the class to be counted and performs the reference counting
// during the life-time of the counted object. When its count goes to zero
// the counted object is destroyed by explicit call to its destructor.
// This class provides overloaded operators *() and ->() to make the
// maniplutation syntax the same as for the normal "dumb" pointers.
// The basic rule for the use of this class is that it is passed always by
// reference and is not constructed by calling new.
// Before the use of this smart pointer object one can test its validity
// using the operator !() or operator bool().
//  if( !smartPtrObj ) { ... } // Problem! We must initialize it first!
//  else               { ... } // OK!
// The code which tries to delete this object won't compile, because it is
// not a pointer it is an object.
// Author:      Radovan Chytracek
// Version:     2.0
// Date:        November 2001
// ----------------------------------------------------------------------
#ifndef _G4REFERENCECOUNTEDHANDLE_H_
#define _G4REFERENCECOUNTEDHANDLE_H_ 1

#include "G4Allocator.hh"
#define RCH_USING_G4ALLOCATOR 1

template <class X> class G4ReferenceCountedHandle
{
public:

  typedef G4ReferenceCountedHandle RCH;
  
private: // with description

  class CountedObject
  {
    
    friend G4ReferenceCountedHandle<X>;
    
  public:

    CountedObject( X* pObj = 0 ) : fCount(0), fRep( pObj ) {
      if( pObj != 0 ) {
        fCount = 1;
      }
    } // Constructor
    
    ~CountedObject() {
      delete fRep;
    } // Destructor
    
    inline void AddRef() {
      ++fCount;
    } // Forward to Counter class
    
    inline void Release() {
      if( --fCount == 0 ) delete this;
    } // Forward to Counter class
    
    typedef G4Allocator<CountedObject> COUNTEDOBJALLOCATOR;
    
    // There is no provision that this class is subclassed.
    // If it is subclassed & new data members are added then the
    // following "new" & "delete" will fail and give errors. 
    inline void* operator new( size_t )	{
      return( (void *)CountedObject::GetAllocator().MallocSingle() );
    } // operator new defined for G4Allocator
    
    inline void operator delete( void *pObj )	{
      CountedObject::GetAllocator().FreeSingle( (CountedObject*)pObj );
    } // operator delete defined for G4Allocator
    
  private:

    unsigned int fCount;
    // Reference counter
    X*           fRep;
    // The counted object
    
  private:

    static COUNTEDOBJALLOCATOR& GetAllocator() {
      return( aCountedObjectAllocator );
    }
    
  private:

    static G4Allocator<CountedObject> aCountedObjectAllocator;
  };
  
public: // with description

  G4ReferenceCountedHandle( X* rep = 0 )	: fObj( 0 ) {
    if( rep != 0 ) {
      fObj = new CountedObject( rep );
    }
  } // Constructor
  
  G4ReferenceCountedHandle( const G4ReferenceCountedHandle& right ) : fObj( right.fObj ) {
    fObj->AddRef();
  } // Copy constructor
  
  ~G4ReferenceCountedHandle()  {
    if( fObj ) fObj->Release();
  } // Destructor
  
  inline G4ReferenceCountedHandle& operator =( const G4ReferenceCountedHandle& right ) {
    if( fObj != right.fObj ) {
      if( fObj )
        fObj->Release();
      this->fObj = right.fObj;
      fObj->AddRef();
    }
    return *this;
  } // Assignment operator by reference
  
  inline G4ReferenceCountedHandle& operator =( X* objPtr ) {
    if( fObj )
      fObj->Release();
    this->fObj = new  CountedObject( objPtr );
    return *this;
  } // Assignment operator by reference
  
  inline unsigned int Count() const {
    return( fObj ? fObj->fCount : 0 );
  } // Forward to Counter class
  
  inline X* operator ->() const {
    return( fObj ? fObj->fRep : 0 );
  } // Operator -> allowing the access to counted object
    // The check for 0-ness is left out for performance reasons, see operator () bellow
    // May be called on initialized smart-pointer only!
  
  inline bool operator !() const {
    return( ( !fObj ) ? true : false );
  } // Validity test operator
  
  inline operator bool() const {
    return( ( fObj ) ? true : false );
  }
  
  inline X* operator ()() const {
    return( fObj ? fObj->fRep : 0 );
  } // Functor operator (for covinience)
  
  typedef  G4Allocator< G4ReferenceCountedHandle<X> > RCHALLOCATOR;
  
  // There is no provision that this class is subclassed.
  // If it is subclassed & new data members are added then the
  // following "new" & "delete" will fail and give errors. 
  inline void* operator new( size_t ) {
    return( (void *)G4ReferenceCountedHandle<X>::GetAllocator().MallocSingle() );
  } // operator new defined for G4Allocator
  
  inline void operator delete( void *pObj )	{
    G4ReferenceCountedHandle<X>::GetAllocator().FreeSingle( (G4ReferenceCountedHandle<X>*)pObj );
  } // operator delete defined for G4Allocator
  
  inline void* operator new( size_t, void *pObj ) {
    return pObj;
  } // This is required by VC++ in order to compile but produces warning about missing operator delete(...)
  // Normally this variant is not needed but when this class is used in the context of STL container
  
private:

  static RCHALLOCATOR& GetAllocator() {
    return aRCHAllocator;
  }
  
private:

  CountedObject*     fObj;
  // Object being the subject to reference counting

private:

  static G4Allocator< G4ReferenceCountedHandle<X> > aRCHAllocator;
};

// In order to save the human's typing and brain the macro is provided for definition of the allocators
#define DEFINE_RCH_ALLOCATOR(CountedType) \
  \
  G4Allocator<G4ReferenceCountedHandle<##CountedType##>::CountedObject> \
  G4ReferenceCountedHandle<##CountedType##>::CountedObject::aCountedObjectAllocator;\
  G4Allocator<G4ReferenceCountedHandle<##CountedType##> > \
  G4ReferenceCountedHandle<##CountedType##>::aRCHAllocator;\
  
#endif // _G4REFERENCECOUNTEDHANDLE_H_

