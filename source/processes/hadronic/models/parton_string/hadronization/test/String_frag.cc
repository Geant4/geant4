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
// $Id: String_frag.cc,v 1.2 2008-12-18 13:02:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//
//             by Gunter Folger, June 1998.
//       class exercising G4FTFModel class.
// ------------------------------------------------------------
#define debout if(false) G4cout
//#define debout if(debprint) G4cout
 
  
#include "G4ExcitedStringDecay.hh"
#include "Randomize.hh"
#include "G4StableIsotopes.hh"

#ifdef Hbook
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/HBookHistogram.h"
#endif
#include "G4Neutron.hh"
#include "G4Proton.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ExcitedMesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ExcitedBaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Quarks.hh"
#include "G4DiQuarks.hh"

#include "G4ParticleTable.hh"

// #include "G4ParticleTable.hh"
// #include "G4LeptonConstructor.hh" 
// #include "G4BaryonConstructor.hh" 
// #include "G4MesonConstructor.hh"
// #include "G4BosonConstructor.hh"
// #include "G4ShortLivedConstructor.hh"
//#include "G4QuarkAnnihilator.hh"

//void ConstructParticle();
G4double Rapidity(const G4double p, const  G4double E);



int main()
{

	G4cout << " Welcome; How many events shal I do ? " << std::flush;
	G4int maxEvents;
	G4cin >> maxEvents;
	
	G4int PDGleft,PDGright;
	G4cout << G4endl << " PDG code for the string ends " << std::flush;
	G4cin >> PDGleft >> PDGright; 
	G4cout <<  " left / right " << PDGleft << " / " << PDGright<<  G4endl;

	
	G4double ptparameter;
	G4cout << G4endl << " sigma pt in GeV " << std::flush;
	G4cin >> ptparameter; 
	G4cout <<  " ptparameter = " << ptparameter <<  G4endl;

	G4double pspin;
	G4cout << G4endl << " pspsin (vector/scalar mesons) " << std::flush;
	G4cin >> pspin; 
	G4cout <<  " pspin = " << pspin <<  G4endl;

	G4cout << G4endl << " String Momentum (GeV) " << std::flush;
        G4double string_momentum;
        G4cin >> string_momentum;
	G4cout <<  " string_momentum = " << string_momentum << " GeV" << G4endl;
         		
	G4cout << G4endl << " String mass (GeV) " << std::flush;
        G4double string_mass;
        G4cin >> string_mass;
	G4cout <<  " string_mass = " << string_mass << " GeV" << G4endl;
         		
	G4cout << G4endl;
	
//----Hbook       
#ifdef Hbook
	HBookFile * hbfile=new HBookFile("Nucleustest2.hbook", 29); 

	HBookTuple t_particles(" secondaries 4Momentum ",1000);
	HBookTuple t_sumP(" Summed 4Momentum ",2000);
	HBookHistogram h_PDGcodes("PDGcodes",20000,-10000.,10000, 1);
	HBookHistogram h_multiplicity("Multiplicity",20, 0.,20., 10);
#endif
//----
	G4VLongitudinalStringDecay * aStringDecayMethod=
//				new G4LundStringFragmentation(ptparameter*GeV);
				new G4LundStringFragmentation();
	aStringDecayMethod->SetVectorMesonProbability(pspin);
	G4VStringFragmentation * stringDecay=new G4ExcitedStringDecay(aStringDecayMethod);

	
	
		
	string_momentum *=GeV;
	string_mass *=GeV;

// ----------- now get all particles of interest ---------
//	ConstructParticle();
   
	G4BosonConstructor Bosons;
	Bosons.ConstructParticle();

	G4MesonConstructor Mesons;
	Mesons.ConstructParticle();

	G4LeptonConstructor Leptons;
	Leptons.ConstructParticle();

	G4BaryonConstructor Baryons;
	Baryons.ConstructParticle();

	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();
	
//----------------------------
	   			  	
	for ( G4int ntimes=0; ntimes< maxEvents; ntimes++) {
	
	   G4bool debprint=ntimes<5;

	
	   G4ParticleMomentum direction(0.,0.,1.);
//	   G4ParticleMomentum direction(G4UniformRand(),G4UniformRand(),G4UniformRand());
	   direction.setMag(1.);
//	   G4ParticleDefinition * proton = G4Proton::Proton();


	   G4Parton *Pleft=new G4Parton(PDGleft);
	   G4Parton *Pright=new G4Parton(PDGright);
	   
	   G4double String_charge=
	   	Pleft->GetDefinition()->GetPDGCharge() +
		Pright->GetDefinition()->GetPDGCharge();
		
	   G4ExcitedString string(Pleft,Pright);
//	   G4ExcitedStringVector strings;
//	   strings.insert(&string);
	   
// very simple creation of e/p

	   G4LorentzVector ptot(G4ThreeVector(0.,0.,string_momentum),
	   		std::sqrt(sqr(string_momentum)+sqr(string_mass)));

	   Pleft->Set4Momentum(0.5*ptot);
	   Pright->Set4Momentum(0.5*ptot);
	   

	   G4KineticTrackVector * result;
	
//	   result=stringDecay->FragmentStrings(&strings);
	   result=aStringDecayMethod->FragmentString(string);

//-------- check charge

	   G4int Charge  = 0;
	   G4int Leptons = 0;
	   G4int Barions = 0;
	   G4double Etot = 0.;
	   G4LorentzVector Psum=0;
	   G4LorentzVector Pparticle;
	   
	   G4cout <<" New event with " << result->size() << "secondaries" << G4endl;
	   for (G4int loop=0;loop< result->size();loop++)
	   {
	        G4ParticleDefinition * pdef=result->operator[](loop)->GetDefinition();
	        Charge += G4int(pdef->GetPDGCharge());
		Leptons+= pdef->GetLeptonNumber();
		Barions+= pdef->GetBaryonNumber();
		Pparticle = result->operator[](loop)->Get4Momentum();
	   	Psum   += Pparticle;
// 	        G4cout << pdef->GetParticleName() << " : "
// 		       << result->operator[](loop)->Get4Momentum()
// 		       << G4endl;

#ifdef Hbook
		h_PDGcodes.accumulate(pdef->GetPDGEncoding());
		h_multiplicity.accumulate(float(result->size()));
		t_particles.column("px",Pparticle.x());
		t_particles.column("py",Pparticle.y());
		t_particles.column("pz",Pparticle.z());
		t_particles.column("e",Pparticle.e());
		
		t_particles.dumpData();
#endif
	   }

// check e/p coservation....
	   G4LorentzVector Pdiff= string.Get4Momentum() - Psum;
	   G4double Ediff = Pdiff.e();
	   G4cout << "energy difference : " << Ediff << G4endl;
	   
#ifdef Hbook
		t_sumP.column("px",Psum.x());
		t_sumP.column("py",Psum.y());
		t_sumP.column("pz",Psum.z());
		t_sumP.column("e",Psum.e());
		
		t_sumP.dumpData();
#endif
	   
//	   G4cout << " total E : "<< Psum.e() << G4endl;
	   
	   if ( std::abs(Charge - String_charge)  > perCent ) 
	                               G4cout << " N Charge Leptons Barions " 
				              << result->size() << " , "
	   				      <<Charge << " , "  
	   				      << Leptons<< " , " 
					      << Barions<< G4endl;

	   std::for_each(result->begin(), result->end(), DeleteKineticTrack());
	   result->clear();
	   delete result;
	}
	
#ifdef Hbook
	hbfile->write();
#endif

// Clean up....
	delete aStringDecayMethod;
	delete  stringDecay;


	return 0;
}

G4double Rapidity(const G4double p, const G4double E)
{
	return 0.5*std::log((E+p)/(E-p));
} 

//**************************************************
// void ConstructParticle()
//   {
//   // In this method, static member functions should be called
//   // for all particles which you want to use.
//   // This ensures that objects of these particle types will be
//   // created in the program. 
// 
//   G4LeptonConstructor Leptons;
//   Leptons.ConstructParticle();
// 
//   G4BosonConstructor Bosons;
//   Bosons.ConstructParticle();
// 
//   G4MesonConstructor Mesons;
//   Mesons.ConstructParticle();
// 
//   G4BaryonConstructor Baryons;
//   Baryons.ConstructParticle();
// 
//   G4ShortLivedConstructor ShortLived;
//   ShortLived.ConstructParticle();
//   }
