// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VStringFragmentation.hh,v 1.1.8.1 1999/12/07 20:51:53 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
#ifndef G4VStringFragmentation_h
#define G4VStringFragmentation_h 1

#include "G4ExcitedString.hh"
//GF#include "G4ReactionProductVector.hh"
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
//      virtual G4ReactionProductVector * FragmentString(const G4ExcitedString &theString)=0;
      virtual G4KineticTrackVector * FragmentString(const G4ExcitedString &theString)=0;

  private:

};

#endif


