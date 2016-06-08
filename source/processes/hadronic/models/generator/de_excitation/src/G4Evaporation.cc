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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Evaporation.cc,v 1.7 2001/08/01 17:05:28 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
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
    theChannels = new G4std::vector<G4VEvaporationChannel*>;
    theChannels->reserve(8);

    theChannels->push_back( new G4NeutronEvaporationChannel() );  // n
    theChannels->push_back( new G4ProtonEvaporationChannel() );   // p
    theChannels->push_back( new G4DeuteronEvaporationChannel() ); // Deuteron
    theChannels->push_back( new G4TritonEvaporationChannel() );   // Triton
    theChannels->push_back( new G4He3EvaporationChannel() );      // He3
    theChannels->push_back( new G4AlphaEvaporationChannel() );    // Alpha

    theChannels->push_back( new G4CompetitiveFission() ); // Fission Channel
    theChannels->push_back( new G4PhotonEvaporation() );  // Photon Channel
}

G4Evaporation::G4Evaporation(const G4Evaporation &right)
{
    G4Exception("G4Evaporation::copy_constructor meant to not be accessable.");
}


G4Evaporation::~G4Evaporation()
{
    if (myOwnChannelsVector) {
	//    theChannels->clearAndDestroy();
	while (!theChannels->empty()) {
	    delete (theChannels->back());
	    theChannels->pop_back();
	}
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
	theResult->push_back(new G4Fragment(theNucleus));
	return theResult;
    }

    // The residual nucleus (after evaporation of each fragment)
    G4Fragment theResidualNucleus = theNucleus;

    // Number of channels
    G4int TotNumberOfChannels = theChannels->size();  
	

    // Starts loop over evaporated particles
    for (;;) {
	// loop over evaporation channels
	G4std::vector<G4VEvaporationChannel*>::iterator i;
	for (i=theChannels->begin(); i != theChannels->end(); i++) {
	    (*i)->Initialize(theResidualNucleus);
	}
	// Work out total decay probability by summing over channels 
	G4double TotalProbability = 0;
	for (i=theChannels->begin(); i != theChannels->end(); i++) {
	    TotalProbability += (*i)->GetEmissionProbability();
	}
	if (TotalProbability <= 0.0) {
	    // Will be no evaporation more
	    // write information about residual nucleus
	    theResult->push_back(new G4Fragment(theResidualNucleus));
	    break; 
	} else {
	    // Selection of evaporation channel, fission or gamma
	    // G4double * EmissionProbChannel = new G4double(TotNumberOfChannels);
	    G4std::vector<G4double> EmissionProbChannel;
	
	    // EmissionProbChannel[0] = theChannels->at(0)->GetEmissionProbability();
	    EmissionProbChannel.push_back(theChannels->front()->GetEmissionProbability()); // index 0
			
	    for (i= (theChannels->begin()+1); i != theChannels->end(); i++) {
		// EmissionProbChannel[i] = EmissionProbChannel[i-1] + 
		// theChannels->at(i)->GetEmissionProbability();
		EmissionProbChannel.push_back(EmissionProbChannel.back() + (*i)->GetEmissionProbability());
	    }

	    G4double shoot = G4UniformRand() * TotalProbability;
	    G4int j;
	    for (j=0; j < TotNumberOfChannels; j++) {
		// if (shoot < EmissionProbChannel[i]) 
		if (shoot < EmissionProbChannel[j]) 
		    break;
	    }
			
	    // delete [] EmissionProbChannel;
	    EmissionProbChannel.clear();
			
	    if( j >= TotNumberOfChannels ) {
		G4Exception( "G4Evaporation::BreakItUp: Can't define emission probability of the channels!" );
	    } else {
		// Perform break-up
		G4FragmentVector * theEvaporationResult = (*theChannels)[j]->BreakUp(theResidualNucleus);

#ifdef debug
		G4cout <<           "---------------------------------------------" 
		       << G4endl;
		G4cout << G4endl << " After evaporate a particle test conservation " 
		       << G4endl;
		CheckConservation(theResidualNucleus,theEvaporationResult);
		G4cout << G4endl 
		       <<           "---------------------------------------------" 
		       << G4endl;
#endif  

		// Check if chosen channel is fission (there are only two EXCITED fragments)
		// or the channel could not evaporate anything
		if ( theEvaporationResult->size() == 1 || 
		     ((*(theEvaporationResult->begin()))->GetExcitationEnergy() > 0.0 && 
		      (*(theEvaporationResult->end()-1))->GetExcitationEnergy() > 0.0) ) {
		    // FISSION 
		    while (theEvaporationResult->size() > 0) {
			theResult->push_back(*(theEvaporationResult->end()-1));
			theEvaporationResult->pop_back();
		    }
		    if (theEvaporationResult->empty()) delete theEvaporationResult;
		    else G4Exception("G4Evaporation.cc: deleting theEvaporationResult but is not empty after fission");
		    break;
		} else {
		    // EVAPORATION
		    while (theEvaporationResult->size() > 1) {
			theResult->push_back(*(theEvaporationResult->begin()));
			theEvaporationResult->erase(theEvaporationResult->begin());
		    }
		    theResidualNucleus = **(theEvaporationResult->begin());
		    delete *(theEvaporationResult->begin());
		    theEvaporationResult->pop_back();
		    if (theEvaporationResult->empty()) delete theEvaporationResult;
		    else G4Exception("G4Evaporation.cc: deleting theEvaporationResult but is not empty after evaporation");
		}
	    }
	}
    }

#ifdef debug
    G4cout << "======== Evaporation Conservation Test ===========" << G4endl 
	   << "==================================================" << G4endl;
    CheckConservation(theNucleus,theResult);
    G4cout << "==================================================" << G4endl;
#endif
    return theResult;
}



#ifdef debug
void G4Evaporation::CheckConservation(const G4Fragment & theInitialState,
				      G4FragmentVector * Result) const
{
    G4double ProductsEnergy =0;
    G4ThreeVector ProductsMomentum;
    G4int ProductsA = 0;
    G4int ProductsZ = 0;
    for (G4FragmentVector::iterator h = Result->begin(); h != Result->end(); h++) {
	G4LorentzVector tmp = (*h)->GetMomentum();
	ProductsEnergy += tmp.e();
	ProductsMomentum += tmp.vect();
	ProductsA += G4int((*h)->GetA());
	ProductsZ += G4int((*h)->GetZ());
    }

    if (ProductsA != theInitialState.GetA()) {
	G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4Evaporation.cc: Barionic Number Conservation test for evaporation fragments" 
	       << G4endl; 
	G4cout << "Initial A = " << theInitialState.GetA() 
	       << "   Fragments A = " << ProductsA << "   Diference --> " 
	       << theInitialState.GetA() - ProductsA << G4endl;
    }
    if (ProductsZ != theInitialState.GetZ()) {
	G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4Evaporation.cc: Charge Conservation test for evaporation fragments" 
	       << G4endl; 
	G4cout << "Initial Z = " << theInitialState.GetZ() 
	       << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	       << theInitialState.GetZ() - ProductsZ << G4endl;
    }
    if (abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
	G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4Evaporation.cc: Energy Conservation test for evaporation fragments" 
	       << G4endl; 
	G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	       << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	       << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
    } 
    if (abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
	abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
	abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
	G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!" << G4endl;
	G4cout << "G4Evaporation.cc: Momentum Conservation test for evaporation fragments" 
	       << G4endl; 
	G4cout << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	       << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	       << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV" << G4endl;
    }
    return;
}
#endif




