// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)

#include "G4Evaporation.hh"

#include "G4NeutronEvaporationChannel.hh"
#include "G4ProtonEvaporationChannel.hh"
#include "G4DeuteronEvaporationChannel.hh"
#include "G4TritonEvaporationChannel.hh"
#include "G4He3EvaporationChannel.hh"
#include "G4AlphaEvaporationChannel.hh"
#include "G4CompetitiveFission.hh"
#include "G4PhotonEvaporation.hh"

G4Evaporation::G4Evaporation() : myOwnChannelsVector(true)
{
	theChannels = new G4RWTPtrOrderedVector<G4VEvaporationChannel>;

	theChannels->insert( new G4NeutronEvaporationChannel() );  // n
	theChannels->insert( new G4ProtonEvaporationChannel() );   // p
	theChannels->insert( new G4DeuteronEvaporationChannel() ); // Deuteron
	theChannels->insert( new G4TritonEvaporationChannel() );   // Triton
	theChannels->insert( new G4He3EvaporationChannel() );      // He3
	theChannels->insert( new G4AlphaEvaporationChannel() );    // Alpha

	theChannels->insert( new G4CompetitiveFission() ); // Fission Channel
	theChannels->insert( new G4PhotonEvaporation() );  // Photon Channel
}

G4Evaporation::G4Evaporation(const G4Evaporation &right)
{
	G4Exception("G4Evaporation::copy_constructor meant to not be accessable.");
}


G4Evaporation::~G4Evaporation()
{
	if (myOwnChannelsVector) {
		theChannels->clearAndDestroy();
		delete theChannels;
	}
}

const G4Evaporation & G4Evaporation::operator=(const G4Evaporation &right)
{
  G4Exception("G4Evaporation::operator= meant to not be accessable.");
  return *this;
}


G4bool G4Evaporation::operator==(const G4Evaporation &right) const
{
  return false;
}

G4bool G4Evaporation::operator!=(const G4Evaporation &right) const
{
  return true;
}


G4FragmentVector * G4Evaporation::BreakItUp(const G4Fragment &theNucleus)
{
	G4FragmentVector * theResult = new G4FragmentVector;

	// CHECK that Excitation Energy != 0
	if (theNucleus.GetExcitationEnergy() <= 0.0) {
		theResult->insert(new G4Fragment(theNucleus));
		return theResult;
	}

	// The residual nucleus (after evaporation of each fragment)
	G4Fragment theResidualNucleus = theNucleus;

	// Number of channels
	G4int TotNumberOfChannels = theChannels->entries();  
	

	// Starts loop over evaporated particles
	for (;;) {
		// loop over evaporation channels
		G4int i;
		for (i=0; i < TotNumberOfChannels; i++) {
			theChannels->at(i)->Initialize(theResidualNucleus);
    	}
		// Work out total decay probability by summing over channels 
		G4double TotalProbability = 0;
		for (i=0; i < TotNumberOfChannels; i++) {
			TotalProbability += theChannels->at(i)->GetEmissionProbability();
		}
		if (TotalProbability <= 0.0) {
			// Will be no evaporation more
			// write information about residual nucleus
			theResult->insert(new G4Fragment(theResidualNucleus));
			break; 
		} else {
			// Selection of evaporation channel, fission or gamma
// 			G4double * EmissionProbChannel = new G4double(TotNumberOfChannels);
			G4RWTValOrderedVector<G4double> EmissionProbChannel;
	
// 			EmissionProbChannel[0] = theChannels->at(0)->GetEmissionProbability();
			EmissionProbChannel.insert(theChannels->at(0)->GetEmissionProbability()); // index 0
			
			for (i=1; i < TotNumberOfChannels; i++) {
// 				EmissionProbChannel[i] = EmissionProbChannel[i-1] + 
// 												 theChannels->at(i)->GetEmissionProbability();
				EmissionProbChannel.insert(EmissionProbChannel(i-1) + 
												 theChannels->at(i)->GetEmissionProbability());
			}

			G4double shoot = G4UniformRand() * TotalProbability;

			for (i=0; i < TotNumberOfChannels; i++) {
// 				if (shoot < EmissionProbChannel[i]) 
				if (shoot < EmissionProbChannel(i)) 
				break;
			}
			
// 			delete [] EmissionProbChannel;
			EmissionProbChannel.clear();
			
			if( i >= TotNumberOfChannels ) {
				G4Exception( "G4Evaporation::BreakItUp: Can't define emission probability of the channels!" );
			} else {
				// Perform break-up
				G4FragmentVector * theEvaporationResult = theChannels->at(i)->BreakUp(theResidualNucleus);

				// Check if chosen channel is fission (there are only two EXCITED fragments)
				// or the channel could not evaporate anything
				if ( theEvaporationResult->entries() == 1 || 
						(theEvaporationResult->first()->GetExcitationEnergy() > 0.0 && 
						theEvaporationResult->last()->GetExcitationEnergy() > 0.0) ) {
					// FISSION 
					while (theEvaporationResult->entries() > 0) theResult->insert(theEvaporationResult->removeFirst());
					theEvaporationResult->clearAndDestroy();
					delete theEvaporationResult;
					break;
				} else {
					// EVAPORATION
					while (theEvaporationResult->entries() > 1) theResult->insert(theEvaporationResult->removeFirst());
					theResidualNucleus = *(theEvaporationResult->at(0));
					theEvaporationResult->clearAndDestroy();
					delete theEvaporationResult;
				}
			}
		}
	}

	return theResult;
}



