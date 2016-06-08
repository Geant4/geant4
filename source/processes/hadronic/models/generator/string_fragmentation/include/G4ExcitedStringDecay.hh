// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedStringDecay.hh,v 1.6 2000/11/10 08:30:07 hpw Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
#ifndef G4ExcitedStringDecay_h
#define G4ExcitedStringDecay_h 1

#include "globals.hh"
#include "G4VStringFragmentation.hh"
#include "G4ExcitedStringVector.hh"
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

      virtual G4KineticTrackVector * FragmentStrings(const G4ExcitedStringVector * theStrings);

  private:
      G4KineticTrackVector * FragmentString(const G4ExcitedString &theString);
      G4bool EnergyAndMomentumCorrector(G4KineticTrackVector* Output, G4LorentzVector& TotalCollisionMom);   
  
      G4VLongitudinalStringDecay * theStringDecay;

};

inline
G4KineticTrackVector *G4ExcitedStringDecay::
FragmentString(const G4ExcitedString &theString)
{
	if ( theStringDecay == NULL ) 
	    theStringDecay=new G4LundStringFragmentation();
	    
	return theStringDecay->FragmentString(theString);
}
	
inline
G4KineticTrackVector *G4ExcitedStringDecay::
FragmentStrings(const G4ExcitedStringVector * theStrings)
{
  G4KineticTrackVector * theResult = new G4KineticTrackVector;

  G4LorentzVector KTsum;
  G4bool NeedEnergyCorrector=false;
  
  for ( G4int astring=0; astring < theStrings->entries(); astring++)
  {
	KTsum+= theStrings->at(astring)->Get4Momentum();
        G4KineticTrackVector * generatedKineticTracks = NULL;
  
	if ( theStrings->at(astring)->IsExcited() )
	{
  	     generatedKineticTracks=FragmentString(*theStrings->at(astring));
	} else {
	     generatedKineticTracks = new G4KineticTrackVector;
	     generatedKineticTracks->insert(theStrings->at(astring)->GetKineticTrack());
	}    

	if (generatedKineticTracks == NULL) 
	{
		G4cerr << "G4VPartonStringModel:No KineticTracks produced" << G4endl;
		continue;
	}
	
	G4LorentzVector KTsum1;
	for ( G4int aTrack=0; aTrack<generatedKineticTracks->entries();aTrack++)
	{
		theResult->insert(generatedKineticTracks->at(aTrack));
		KTsum1+= (*generatedKineticTracks)[aTrack]->Get4Momentum();
	}
	
	
	if  ( abs((KTsum1.e()-theStrings->at(astring)->Get4Momentum().e()) / KTsum1.e()) > perMillion ) 
	{
	   NeedEnergyCorrector=true;
 	}
//      clean up
	delete generatedKineticTracks;
  }
//  G4cout << "String total energy and 4 momentum" <<KTsum.t()<<" "<< KTsum << endl;
  
  if ( NeedEnergyCorrector ) EnergyAndMomentumCorrector(theResult, KTsum);
    
  return theResult;
}

#endif


