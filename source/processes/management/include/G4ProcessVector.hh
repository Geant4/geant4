// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessVector.hh,v 1.5 1999-11-15 10:39:47 gunter Exp $
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
//
// Class Description
//  This class is a container for pointers to physics process objects. Its 
//  functionality is derived from G4RWTPtrOrderedVector<T>.
// ------------------------------------------------------------

#ifndef G4ProcessVector_h
#define G4ProcessVector_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4rw/tpordvec.h"

#include "G4VProcess.hh"

class G4ProcessVector : public G4RWTPtrOrderedVector<G4VProcess>
{
  //  Is a container for pointers to physics process objects. Its 
  //  functionality is derived from G4RWTPtrOrderedVector<T>.

  public:

      G4ProcessVector(size_t capac=G4RWDEFAULT_CAPACITY)
        : G4RWTPtrOrderedVector<G4VProcess>(capac) {;}
      //  Constructor.

      virtual ~G4ProcessVector() {;}
      //  Destructor.

      G4bool operator==(const G4ProcessVector &right) const
      {
        return (this == (G4ProcessVector *) &right);
      }

};
  
#endif
