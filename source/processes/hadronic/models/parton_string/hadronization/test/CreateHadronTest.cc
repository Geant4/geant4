// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: CreateHadronTest.cc,v 1.1 2003-10-07 11:27:10 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//
//             by Gunter Folger, June 1998.
//       class exercising small part of G4VLongitudinalStringDecay.
// ------------------------------------------------------------
	   
#include "G4PionPlus.hh"
#include "G4ReactionProductVector.hh"
#include "G4ExcitedStringDecay.hh"
#include "Randomize.hh"
#include "G4StableIsotopes.hh"


#include "G4ParticleTable.hh"
#include "G4LeptonConstructor.hh" 
#include "G4BaryonConstructor.hh" 
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
//#include "G4QuarkAnnihilator.hh"


void ConstructParticle();

int main()
{

	G4cout << " Welcome; How many events shal I do ? " ;
	G4cout.flush();
	G4int maxEvents;
	G4cin >> maxEvents;
	
	G4cout << " Quark 1 + 2 PDG codes " ;
	G4cout.flush();
	G4int q1,q2;
	G4cin >> q1 >> q2;
	
	G4cout << " Spin L (2J+1 ), -1 random " ;
	G4cout.flush();
	G4int L;
	G4cin >> L;
	G4cout << endl;

       
	ConstructParticle();
	
	G4VLongitudinalStringDecay * aStringDecay=
				new G4LundStringFragmentation();

	for (G4int i=0; i<=maxEvents;  ++i)
	{
	    G4ParticleDefinition * particle=
	              aStringDecay->CreateHadron(q1,q2, L>= 0, L);
	    G4cout << particle->GetPDGEncoding() << " , "
		   << particle->GetParticleName();
	    G4double quarkcharge=
		G4ParticleTable::GetParticleTable()->FindParticle(q1)->GetPDGCharge() +
	    	G4ParticleTable::GetParticleTable()->FindParticle(q2)->GetPDGCharge();
	    G4double charge_balance=particle->GetPDGCharge() - quarkcharge;
	    if (abs(charge_balance) > 0.01 ) G4cout << " incorrect charge: "
		    				    << charge_balance;
	    G4cout   << G4endl;
	    
	    
	}
}
//**************************************************
void ConstructParticle()
  {
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();

  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4BaryonConstructor Baryons;
  Baryons.ConstructParticle();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  }
