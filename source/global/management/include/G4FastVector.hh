// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastVector.hh,v 1.1 1999-01-07 16:09:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// ------------------------------------------------------------

#ifndef G4FastVector_h
#define G4FastVector_h 1

#include "globals.hh"
#include "G4ios.hh"


template <class Type, G4int N>
class G4FastVector 
{
  //  Template class defining a vector of pointers,
  //  not performing boundary checking.

  public:

      G4FastVector() { ptr = &theArray[0]; }

      ~G4FastVector()
      {
        if (ptr != &theArray[0]) delete [] ptr;
      }

      inline Type* operator[](G4int anIndex) const
      //  Access operator to the array.
      {
        return ptr[anIndex];
      }

      void Initialize(G4int items)
      //  Normally the pointer ptr points to the stack-array
      //  theArray; only when the number of items is greater
      //  than N, memory is allocated dynamically.
      {
        if (ptr != &theArray[0])
           delete [] ptr;
        if (items > N)
           ptr = new Type*[items];
        else
           ptr = &theArray[0];
      } 

      inline void SetElement(G4int anIndex, Type *anElement)
      //  To insert an element at the given position inside
      //  the vector.
      {
        ptr[anIndex] = anElement;
      }

  private:

      Type *theArray[N];
      Type **ptr;
};

#endif

