// Class RCObjectHandle
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
// Version:     1.0
// Date:        June 2001
// ----------------------------------------------------------------------
#ifndef RC_OBJECT_HANDLE_H
#define RC_OBJECT_HANDLE_H 1

template <class X>
class RCObjectHandle
{
public:
  typedef RCObjectHandle RCH;

private: // with description
    
  class CountedObject {
  public:
    inline CountedObject( X* pObj = 0 )
    : fCount(1), fRep( pObj ) {
    } // Constructor

    inline ~CountedObject() {
      delete fRep;
      fRep = 0;
    } // Destructor
    
    inline X* GetObject() const {
      return fRep;
    }
    
    inline void SetObject( const X* object ) {
      fRep = object;
    }
    
    inline unsigned int AddRef() {
       return( ++fCount );
    } // Forward to Counter class

    inline unsigned int Release() {
       return( ( fCount > 0 ) ? --fCount : 0 );
    } // Forward to Counter class

    inline unsigned int Count() const {
       return( fCount );
    } // Forward to Counter class
    
  private:
    unsigned int fCount;
    // Reference counter
    X*           fRep;
    // The counted object
  };

public: // with description
  RCObjectHandle( X* rep = 0 )
  : fObj( 0 ) {
    if( rep != 0 ) {
	    fObj = new RCH::CountedObject( rep );
	  }
  } // Constructor

  ~RCObjectHandle() {
    if( 0 == Release() )    {
      if( fObj != 0 )     {
        delete fObj;
	      fObj = 0;
      }
    }
  } // Destructor

  RCObjectHandle( const RCObjectHandle& right )
  : fObj( 0 ) {
    fObj = right.GetCountedObject();
    AddRef();
  } // Copy constructor
  
  template <class T>
  RCObjectHandle( const RCObjectHandle<T>& right )
  : fObj( 0 ) {
    fObj = right.GetCountedObject();
    AddRef();
  } // Copy constructor

  RCObjectHandle& operator =( const RCObjectHandle& right ) {
    if( this->fObj != right.GetCountedObject() )              {
      if( 0 == this->Release() )          {
        if( this->fObj != 0 )           {
          delete this->fObj;
        }
      }
      this->fObj = right.GetCountedObject();
      this->AddRef();
    }
    return *this;
  } // Assignment operator by reference
  
  template <class T>
  RCObjectHandle& operator=( const RCObjectHandle<T>& right ) {
    if( this->fObj != right.GetCountedObject() )              {
      if( 0 == this->Release() )          {
        if( this->fObj != 0 )           {
          delete this->fObj;
        }
      }
      this->fObj = right.GetCountedObject();
      this->AddRef();
    }
    return *this;
  } // Assignment operator by reference

  RCObjectHandle& operator =( const RCObjectHandle* right ) {
    if( this->fObj != right->GetCountedObject() )               {
      if( 0 == this->Release() )          {
        if( this->fObj != 0 )           {
          delete this->fObj;
        }
      }
      this->fObj = right->GetCountedObject();
      this->AddRef();
    }
    return *this;
  } // Assignment operator by pointer, should not be ever called 

  RCObjectHandle& operator =( X* pRefObj ) {
    if( this->fObj == 0 || pRefObj != this->fObj->GetObject() ) {
      if( 0 == Release() )          {
        if( fObj != 0 )           {
          delete fObj;
        }
      }
      fObj = new RCH::CountedObject( pRefObj );
    }
    return *this;
  } // Assignment operator by pointer to the counted class object
  
  RCH::CountedObject* GetCountedObject() const {
    return fObj;
  }

  inline unsigned int AddRef() {
     return( ( fObj != 0 ) ? fObj->AddRef() : 0  );
  } // Forward to Counter class

  inline unsigned int Release() {
     return( (fObj != 0 ) ? fObj->Release() : 0 );
  } // Forward to Counter class

  inline unsigned int Count() const {
     return( ( fObj != 0 ) ? fObj->Count() : 0 );
  } // Forward to Counter class

  X* operator ->() const {
     return( ( fObj != 0 ) ? fObj->GetObject() : 0 );
  } // Operator -> allowing the access to counted object
    // The check for 0-ness is left out for performance reasons, see operator () bellow
    // May be called on initialized smart-pointer only!

  bool operator !() const {
     return( ( fObj == 0 || fObj->GetObject() == 0 ) ? true : false );
  } // Validity test operator
  
  operator bool() const {
    return( ( fObj != 0 && fObj->GetObject() != 0 ) ? true : false );
  }

  X& operator *() const {
     return( ( fObj != 0 ) ? *(fObj->GetObject()) : 0 );
  } // Dereference operator to make the feeling of dereferencing a pointer to
    // the counted object
    // May be called on initialized smart-pointer only!

  X* operator ()() const {
     return( ( fObj != 0 ) ? fObj->GetObject() : 0 );
  } // Functor operator (for covinience)


private:
  RCH::CountedObject*     fObj;
 // Object being the subject to reference counting
};

#endif // RC_OBJECT_HANDLE_H

