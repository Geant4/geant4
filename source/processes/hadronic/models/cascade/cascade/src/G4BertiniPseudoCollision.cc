#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4Nucleus.hh"

#include "G4BertiniPseudoCollision.hh"

G4BertiniPseudoCollision::G4BertiniPseudoCollision() {
  ;
}

G4BertiniPseudoCollision::~G4BertiniPseudoCollision(){
  ;
}

/*
	// pseudo collision 1
	++pdevts[mat];
	++pdevst;
	if (tip[no] > 1.0) {
	  if (ec[no] > ehipi) ++hpdevt[mat];
	} else {
	  if (ec[no] > ehin) ++hpdevt[mat];
	}
	lelem = 0;
	nopart = -1;

// pseudo collision 2
	    ++phevts[mat];
	    if (ihie == 1) ++hpevth[mat];
	    iphev  = 1;
	    lelem  = 0;
	    nopart = -1;
*/


