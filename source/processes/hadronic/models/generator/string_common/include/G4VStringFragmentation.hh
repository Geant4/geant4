// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VStringFragmentation.hh,v 1.3 2000/08/02 08:13:05 hpw Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
#ifndef G4VStringFragmentation_h
#define G4VStringFragmentation_h 1

#include "G4ExcitedStringVector.hh"

class G4KineticTrackVector;

class G4VStringFragmentation 
{
  public:
      G4VStringFragmentation();
      ~G4VStringFragmentation();

  private:
      G4VStringFragmentation(const G4VStringFragmentation &right);
      const G4VStringFragmentation & operator=(const G4VStringFragmentation &right);
      int operator==(const G4VStringFragmentation &right) const;
      int operator!=(const G4VStringFragmentation &right) const;

  public:
      virtual G4KineticTrackVector * FragmentStrings(const G4ExcitedStringVector * theStrings)=0;

  private:

};

#endif


