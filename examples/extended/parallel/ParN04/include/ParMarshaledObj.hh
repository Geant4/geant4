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
// $Id: ParMarshaledObj.hh,v 1.1 2002-03-05 15:22:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------

#ifndef ParMarshaledObj_h
#define ParMarshaledObj_h

#include <stdlib.h>
#include "g4std/iostream"

// Geant4 types
#include "globals.hh"
#include "G4Event.hh"

// A MarshaledObj is conceptually a
//                     dynamically resizable buffer, size, and cursor.
//   Data structure of buffer:  <field1, field2, ...>
//   Note that buffer has no internal construct to determine size of buffer
//   It is the responsibility of the user to know types of field1, ...
// There are two types of MarshaledObj, depending on whether
//   IsUnmarshaling() is true or false.
// I use malloc/free because I can also use realloc to extend the buffer.

class MarshaledObj {
  private:
    // Make sure all marshaled objects are word aligned.
    static const int WORD_SIZE = sizeof(long);
    static inline int ROUND_UP( int x )
    { return (((x)+(WORD_SIZE-1)) / WORD_SIZE) * WORD_SIZE; }

  public:
    inline MarshaledObj(int total_size)
      :extent(ROUND_UP(total_size)),size(0)
      { if (extent > 0)
          buffer = (char *)malloc(extent);
        else
          buffer = NULL;
        cursor = buffer + size;
      }
    inline MarshaledObj()
      :extent(128),size(0)
      { buffer = (char *)malloc(extent);
        cursor = buffer + size;
      }
    // this is for an arbitrary buffer; we specify the size explicitly
    // if buf_len == 0, then we'll be unmarshaling
    inline MarshaledObj(void *buf, int buf_len = 0)
      :extent(ROUND_UP(buf_len+sizeof(int))),
	size(ROUND_UP(buf_len+sizeof(buf_len)))
      { if ( buf_len <= 0 ) {
          extent = size = 0; // IMPLIES:  IsUnmarshaling() => true
          buffer = (char *)buf;
          cursor = buffer;
        } else {
          buffer = (char *)malloc(extent);
          *(int *)buffer = buf_len;
          memcpy( buffer+sizeof(int), (char *)buf, extent );
          cursor = buffer + size;
        }
      }
    inline ~MarshaledObj()
    { if ( ! IsUnmarshaling() )
        free(buffer);
    }

    inline bool IsUnmarshaling() { return (extent <= 0); }

  private:
    // Don't use copy constructor
    const MarshaledObj& operator=(const MarshaledObj& right);

  protected:
    char *buffer; /* First entry is size of buffer */
  private:
    int extent;   /* extent is number of bytes allocated, extent >= size */
    int size;     /* Initially, buffer has entry for size */
    char *cursor;

    inline void ResizeBuffer( size_t new_size )
      { int displacement = cursor - buffer;
        extent = new_size;
        buffer = (char *)realloc( buffer, extent );
        cursor = buffer + displacement;
      }

  public:
    inline int BufferSize()
      { return cursor - buffer; }
    inline void *GetBuffer() {
      return buffer;
    }
    // inline void SetNewBuffer(void *buf) {
    //   free(buffer);
    //   buffer = (char *)buf;
    //   size = 0; /* unknown size */
    //   cursor = buffer;
    // }

#define MarshalDefinition(TYPE) \
    inline void Marshal(TYPE a##TYPE) { \
      if (IsUnmarshaling()) \
        throw "Tried to marshal in object marked isUnmarshaling = true"; \
      size += ROUND_UP(sizeof(a##TYPE)); \
      while ( size > extent ) \
        ResizeBuffer( 2 * extent ); \
      *(TYPE *)cursor = a##TYPE; \
      cursor = buffer + size; \
    } \
    inline void Unmarshal(TYPE & a##TYPE) { \
      a##TYPE = *(TYPE *)cursor; \
      cursor += ROUND_UP(sizeof(a##TYPE)); \
    } \
    inline TYPE UnmarshalAs##TYPE() { \
      TYPE a##TYPE = *(TYPE *)cursor; \
      cursor += ROUND_UP(sizeof(a##TYPE)); \
      return a##TYPE; \
    }

// This works for any inline objects (no embedded pointers)
MarshalDefinition(G4int)
MarshalDefinition(G4double)
MarshalDefinition(G4ThreeVector)
// MarshalDefinition(int)
// MarshalDefinition(double)


// These cannot follow the pattern, since the data is not inline.
    inline void MarshalG4String(G4String& aG4String) {
      if (IsUnmarshaling())
        throw "Tried to marshal in object marked isUnmarshaling = true";
      size += ROUND_UP( aG4String.length() + 1 ); // 1 for final NULL char
      while ( size > extent )
        ResizeBuffer( 2 * extent );
      strcpy( (char *)cursor, aG4String );
      cursor = buffer + size;
    }
    inline void Marshal(G4String& aG4String)
    { MarshalG4String( aG4String ); }
    inline void Unmarshal(G4String& aG4String) {
      aG4String = (char *)cursor;
      cursor += ROUND_UP( aG4String.length() + 1 ); // 1 for final NULL char
    }
    // INEFFICIENT? (copies return value)
    inline G4String UnmarshalAsG4String() {
      G4String aG4String = (char *)cursor;
      cursor += ROUND_UP( aG4String.length() + 1 ); // 1 for final NULL char
      return aG4String;
    }

// USAGE:  MarshalObjPtr( &aG4obj, sizeof(aG4obj) );
// USAGE:  MarshaledObj tmp;
//         G4Obj *aObjPtr = tmp.UnmarshalObjPtr();
//          ...
//         delete aObjPtr;
    inline void MarshalObjPtr(const void * aObjPtr, int obj_size)
    {  
      if (IsUnmarshaling())
        throw "Tried to marshal in object marked isUnmarshaling = true";
      size += ROUND_UP( obj_size + sizeof(int) );
      while ( size > extent )
        ResizeBuffer( 2 * extent );
      *(int *)cursor = obj_size;
      memcpy( cursor+sizeof(int), aObjPtr, obj_size );
      cursor = buffer + size;
    }
    inline void Marshal(const void * aObjPtr, int obj_size)
    { MarshalObjPtr(aObjPtr,obj_size); }
    inline void Unmarshal(void * & aObjPtr)
    {
      int obj_size = *(int *)cursor;
      aObjPtr = malloc( obj_size );
      memcpy( aObjPtr, cursor+sizeof(int), obj_size );
      cursor += ROUND_UP( obj_size + sizeof(int) );
    }
    inline void * UnmarshalAsObjPtr()
    { 
      void *aObjPtr;
      Unmarshal(aObjPtr);
      return aObjPtr;
    }
    inline void * UnmarshalAsSharedObjPtr()
    { 
      int obj_size = *(int *)cursor;
      void *aObjPtr;
      aObjPtr = cursor+sizeof(int);
      cursor += ROUND_UP( obj_size + sizeof(int) );
      return aObjPtr;
    }
    // UnmarshalAsSharedObjPtr:  return value points into buffer of curr. obj.

    inline void MarshalMarshaledObj(MarshaledObj& aMarshaledObj) {
      if (IsUnmarshaling())
        throw "Tried to marshal in object marked isUnmarshaling = true";
      MarshalObjPtr( aMarshaledObj.GetBuffer(), aMarshaledObj.BufferSize() );
    }
    // INEFFICIENT? (copies return value)
    inline MarshaledObj UnmarshalAsMarshaledObj() {
      MarshaledObj aMarshaledObj( ((char *)cursor)+sizeof(int), *(int *)cursor );
      cursor += ROUND_UP( aMarshaledObj.BufferSize()+sizeof(int) );
      return aMarshaledObj;
    }
};

class MarshaledHCofThisEvent : public MarshaledObj
{
  public:
    inline MarshaledHCofThisEvent(int size) : MarshaledObj(size) { }
    // Trivial constructor.  MarshaledHCofThisEvent(0) creates nothing.
    MarshaledHCofThisEvent();
    // Constructor copies hits from HCE of master in this marshaled object
    MarshaledHCofThisEvent(G4Event *anEvent);
    // Constructor copies hits from HCE of anEvent in this marshaled object
    inline MarshaledHCofThisEvent(void *buf)
    : MarshaledObj(buf)  // Object will be set for IsUnmarshaling
    { }

    void UnmarshalSlaveHCofThisEvent();
    // Unmarshal copies hits from slave in HCE of master
    // We inherit destructor of ParMarshaledObj
    // It does so by calling UnmarshalSlaveHitsCollection for dispatching
    void UnmarshalSlaveHitsCollection(G4String & HCname, G4String & SDname);

    void MarshalHitsCollection(G4VHitsCollection * aBaseHC, G4String HCname);
    // To be defined in ParExN02MarshaledHits.cc or other
    //   application-specific file.

  private:
    // Don't use copy constructor
    const MarshaledHCofThisEvent& operator=(const MarshaledHCofThisEvent& right);
    void MarshalHitsCollection(const G4Event *anEvent);
};

/*
main() {
  MarshaledObj tmp = MarshaledObj();
  tmp.Marshal(5);
  tmp.Marshal(5.3);
  tmp.Marshal(7);

  MarshaledObj tmp2;
  tmp2.MarshalMarshaledObj(tmp);
  cout << tmp2.UnmarshalAsG4int() << "\n";
  cout << tmp2.UnmarshalAsG4double() << "\n";
  cout << tmp2.UnmarshalAsG4int() << "\n";
}
*/

#endif
