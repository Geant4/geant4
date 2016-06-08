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
//
// $Id: G4ExcitedStringDecay.hh,v 1.10 2002/12/12 19:17:56 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
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
  
  for ( unsigned int astring=0; astring < theStrings->size(); astring++)
  {
	KTsum+= theStrings->operator[](astring)->Get4Momentum();
        G4KineticTrackVector * generatedKineticTracks = NULL;
  
	if ( theStrings->operator[](astring)->IsExcited() )
	{
  	     generatedKineticTracks=FragmentString(*theStrings->operator[](astring));
	} else {
	     generatedKineticTracks = new G4KineticTrackVector;
	     generatedKineticTracks->push_back(theStrings->operator[](astring)->GetKineticTrack());
	}    

	if (generatedKineticTracks == NULL) 
	{
		G4cerr << "G4VPartonStringModel:No KineticTracks produced" << G4endl;
		continue;
	}
	
	G4LorentzVector KTsum1;
	for ( unsigned int aTrack=0; aTrack<generatedKineticTracks->size();aTrack++)
	{
		theResult->push_back(generatedKineticTracks->operator[](aTrack));
		KTsum1+= (*generatedKineticTracks)[aTrack]->Get4Momentum();
	}
	
	
	if  ( abs((KTsum1.e()-theStrings->operator[](astring)->Get4Momentum().e()) / KTsum1.e()) > perMillion ) 
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


