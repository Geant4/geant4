// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReferenceCountedHandle.hh,v 1.2 2001-03-14 07:20:49 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4RereferenceCountedHandle
//
// Class description:
//
// A class to provide reference counting mechanism for a class which
// is designed without this support. It is a templated class, a smart pointer,
// wrapping the class to be counted and performs the reference counting
// during the life-time of the counted object. When its count goes to zero
// the counted object is destroyed by explicit call to its destructor.
// This class provides overloaded operators *() and ->() to make the
// maniplutation syntax the same as for the normal "dumb" pointers.
// The basic rule for the use of this class is that it is passed always by
// reference and is not constructed by calling new.
//
// Before the use of this smart pointer object one can test its validity
// using the operator !() or operator bool().
//  if( !smartPtrObj ) { ... } // Problem! We must initialize it first!
//  else               { ... } // OK!
//
// The code which tries to delete this object won't compile, because it is
// not a pointer it is an object.
//
// Author:      Radovan Chytracek
// Version:     1.0
// Date:        February 2001
// ----------------------------------------------------------------------
#ifndef _G4REFERENCECOUNTEDHANDLE_H_
#define _G4REFERENCECOUNTEDHANDLE_H_ 1

template <class X>
   class G4ReferenceCountedHandle
   {
   public: // with description
      G4ReferenceCountedHandle( X* rep = 0 )
      : fObj( new G4ReferenceCountedHandle::CountedObject(rep) ) {
      } // Constructor
   
      ~G4ReferenceCountedHandle()       {
         if( 0 == Release() )       {
            if( fObj->fRep != 0 ) {
               delete fObj;
            }
         }
      } // Destructor
   
      G4ReferenceCountedHandle( const G4ReferenceCountedHandle& right )
      : fObj( new G4ReferenceCountedHandle::CountedObject() ) {
         fObj = right.fObj;
         AddRef();
      } // Copy constructor
      
      template <class T>
      G4ReferenceCountedHandle( const G4ReferenceCountedHandle<T>& right )
      : fObj( new G4ReferenceCountedHandle::CountedObject() ) {
         fObj = right.fObj;
         AddRef();
      } // Copy constructor
   
      G4ReferenceCountedHandle& operator =( const G4ReferenceCountedHandle& right ) {
         if( &right != this )            {
            if( 0 == Release() ) {
               if( fObj->fRep != 0 ) {
                  delete fObj;
               }
            }
            fObj = right.fObj;
            AddRef();
         }
         return *this;
      } // Assignment operator by reference
      
      template <class T>
      G4ReferenceCountedHandle& operator =( const G4ReferenceCountedHandle<T>& right ) {
         if( &right != this )            {
            if( 0 == Release() ) {
               if( fObj->fRep != 0 ) {
                  delete fObj;
               }
            }
            fObj = right.fObj;
            AddRef();
         }
         return *this;
      } // Assignment operator by reference
   
      G4ReferenceCountedHandle& operator =( const G4ReferenceCountedHandle* right ) {
         if( this != right )           {
            if( 0 == Release() )      {
               if( fObj->fRep != 0 ) {
                  delete fObj;
               }
            }
            fObj = right->fObj;
            AddRef();
         }
         return *this;
      } // Assignment operator by pointer
   
      G4ReferenceCountedHandle& operator =( const X* pRefObj ) {
         if( 0 == Release() )      {
            if( fObj->fRep != 0 ) {
               delete fObj;
            }
         }
         fObj = new G4ReferenceCountedHandle::CountedObject( pRefObj );
         return *this;
      } // Assignment operator by pointer to the counted class object
   
      inline unsigned int AddRef() {
         return ++fObj->fCount;
      } // Forward to Counter class
   
      inline unsigned int Release() {
         return --fObj->fCount;
      } // Forward to Counter class
   
      inline unsigned int Count() const {
         return fObj->fCount;
      } // Forward to Counter class
   
      X* operator ->() const {
         return fObj->fRep;
      } // Operator -> allowing the access to counted object
   
      bool operator !() const {
         return fObj->fRep == 0 ? true : false;
      } // Validity test operator
      
      operator bool() const {
        return fObj->fRep != 0 ? true : false;
      }
   
      X& operator *() const {
         return *fObj->fRep;
      } // Dereference operator to make the feeling of dereferencing a pointer to
        // the counted object
   
      X* operator ()() const {
         return fObj->fRep;
      } // Functor operator (for covinience)
   
   private: // with description
      struct CountedObject {
         unsigned int fCount;
         // Reference counter
         X*           fRep;
         // The counted object
      
         CountedObject::CountedObject( X* pObj = 0 )
         : fCount(1), fRep( pObj ) {
         } // Constructor
         
         CountedObject::~CountedObject() {
            delete fRep;
         } // Destructor
      };
   
   private:
      G4ReferenceCountedHandle::CountedObject*     fObj;
     // Object being the subject to reference counting
   };

#endif // _G4REFERENCECOUNTEDHANDLE_H_
