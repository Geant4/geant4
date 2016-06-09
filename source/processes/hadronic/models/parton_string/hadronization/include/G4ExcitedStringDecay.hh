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
// $Id: G4ExcitedStringDecay.hh,v 1.7 2007/05/03 22:06:17 gunter Exp $
// GEANT4 tag $Name: geant4-08-03 $
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
  G4LorentzVector KTsecondaries;
  G4bool NeedEnergyCorrector=false;
  
  for ( unsigned int astring=0; astring < theStrings->size(); astring++)
  {
	KTsum+= theStrings->operator[](astring)->Get4Momentum();
	if( !(KTsum.e()<1) && !(KTsum.e()>-1) )
	{
          throw G4HadronicException(__FILE__, __LINE__, 
	                           "G4ExcitedStringDecay::FragmentStrings received nan string...");
	}
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

	KTsecondaries+=KTsum1;
	
	if  ( KTsum1.e() > 0 && std::abs((KTsum1.e()-theStrings->operator[](astring)->Get4Momentum().e()) / KTsum1.e()) > perMillion ) 
	{
//--debug--           G4cout << "String secondaries(" <<generatedKineticTracks->size()<< ")  momentum: " 
//--debug--	          << theStrings->operator[](astring)->Get4Momentum() << " " << KTsum1 << G4endl;
	   NeedEnergyCorrector=true;
 	}
//      clean up
	delete generatedKineticTracks;
  }
//--DEBUG  G4cout << "Strings/secs total  4 momentum " << KTsum << " " <<KTsecondaries << G4endl;

  G4bool success=true;
  if ( NeedEnergyCorrector ) success=EnergyAndMomentumCorrector(theResult, KTsum);


#ifdef debug_ExcitedStringDecay
  G4LorentzVector  KTsum1=0;
  for ( unsigned int aTrack=0; aTrack<theResult->size();aTrack++)
  {
      G4cout << " corrected tracks .. " << (*theResult)[aTrack]->GetDefinition()->GetParticleName()
      <<"  " << (*theResult)[aTrack]->Get4Momentum() << G4endl;
      KTsum1+= (*theResult)[aTrack]->Get4Momentum();
  }
  G4cout << "Needcorrector/success " << NeedEnergyCorrector << "/" << success << ", Corrected total  4 momentum " << KTsum1  << G4endl;
  if ( ! success ) G4cout << "failed to correct E/p" << G4endl;  
#endif

  return theResult;
}

#endif


