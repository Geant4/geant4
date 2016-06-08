// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedStringDecay.hh,v 1.1.10.1 1999/12/07 20:51:55 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#ifndef G4ExcitedStringDecay_h
#define G4ExcitedStringDecay_h 1

#include "G4VStringFragmentation.hh"
#include "G4ExcitedString.hh"
#include "G4KineticTrackVector.hh"
#include "G4LundStringFragmentation.hh"

class G4ExcitedStringDecay: public G4VStringFragmentation 
{
  public:
      G4ExcitedStringDecay();
      G4ExcitedStringDecay(G4VLongitudinalStringDecay * aStringDecay);
      ~G4ExcitedStringDecay();

  private:
      G4ExcitedStringDecay(const G4ExcitedStringDecay &right);
      const G4ExcitedStringDecay & operator=(const G4ExcitedStringDecay &right);
      int operator==(const G4ExcitedStringDecay &right) const;
      int operator!=(const G4ExcitedStringDecay &right) const;

  public:
//      virtual G4ReactionProductVector * FragmentString(const G4ExcitedString &theString)=0;
      virtual G4KineticTrackVector * FragmentString(const G4ExcitedString &theString);

  private:
  
      G4VLongitudinalStringDecay * theStringDecay;

};

inline
G4KineticTrackVector *G4ExcitedStringDecay::FragmentString(const G4ExcitedString &theString)
{
	if ( theStringDecay == NULL ) 
	    theStringDecay=new G4LundStringFragmentation();
	    
	return theStringDecay->FragmentString(theString);
}
	
#endif


