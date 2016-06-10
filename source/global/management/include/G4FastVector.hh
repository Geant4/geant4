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
//
// $Id: G4FastVector.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// ------------------------------------------------------------

#ifndef G4FastVector_h
#define G4FastVector_h 1

#include "globals.hh"

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

