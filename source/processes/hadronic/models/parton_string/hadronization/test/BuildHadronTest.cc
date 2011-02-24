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
//       class exercising G4HadronBuilder.cc
//             created: Gunter Folger, September 2001.
// ------------------------------------------------------------
	   
//#include "G4PionPlus.hh"
//#include "G4ReactionProductVector.hh"
//#include "G4ExcitedStringDecay.hh"
#include "Randomize.hh"
//#include "G4StableIsotopes.hh"
#include "G4HadronBuilder.hh"
#include "G4Parton.hh"


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

	G4double barionMix= 0.5;
	G4double mesonMix= 0.5;
	
	std::vector<G4double> scalarMesonMix(6);
	scalarMesonMix[0]=0.5;
	scalarMesonMix[1]=0.25;
	scalarMesonMix[2]=0.5;
	scalarMesonMix[3]=0.25;
	scalarMesonMix[4]=1.0;
	scalarMesonMix[5]=0.5;  
	
	std::vector<G4double> vectorMesonMix(6);
	vectorMesonMix[0]=0.5;
	vectorMesonMix[1]=0.;
	vectorMesonMix[2]=0.5;
	vectorMesonMix[3]=0.;
	vectorMesonMix[4]=1.0;
	vectorMesonMix[5]=1.0;
	
	for (G4int im=0; im<6; im++) {
		G4cout << "scalar/vector Meson mixing : " << 
				scalarMesonMix[im] << " / " <<
				vectorMesonMix[im] << G4endl;
	}

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
	G4int Spin;
	G4cin >> Spin;
	G4cout << endl;

       
	ConstructParticle();
	
	G4HadronBuilder * hadron = new G4HadronBuilder(mesonMix,barionMix,
						scalarMesonMix,vectorMesonMix);

	G4Parton *parton1=new G4Parton(q1);
	G4Parton *parton2=new G4Parton(q2);
	G4ParticleDefinition * particle;
	
	for (G4int i=0; i<=maxEvents;  ++i)
	{
	    
	    if ( Spin == -1 ) {
		particle=hadron->Build(parton1,parton2);
	    } else if (Spin == 1 || Spin == 2) {
		particle=hadron->BuildLowSpin(parton1,parton2);
	    } else if (Spin == 3 || Spin == 4){
		particle=hadron->BuildHighSpin(parton1,parton2);
	    } else {
		G4cout << " invalid Spin given ---> STOP" << endl;
		exit(1);
	    } 
		    
	    G4cout << particle->GetPDGEncoding() << " , "
		   << particle->GetParticleName();
	    G4double quarkcharge=
		G4ParticleTable::GetParticleTable()->FindParticle(q1)->GetPDGCharge() +
	    	G4ParticleTable::GetParticleTable()->FindParticle(q2)->GetPDGCharge();
	    G4double charge_balance=particle->GetPDGCharge() - quarkcharge;
	    if (std::abs(charge_balance) > 0.01 ) G4cout << " incorrect charge: "
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
