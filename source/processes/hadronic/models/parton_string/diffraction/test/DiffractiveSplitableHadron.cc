// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: DiffractiveSplitableHadron.cc,v 1.1 2003-10-08 13:48:52 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//
//             by Gunter Folger, June 1998.
//       class exercising G4DiffractiveSplitableHadron.
// ------------------------------------------------------------
#include "G4DiffractiveSplitableHadron.hh"

#include "G4ParticleTable.hh"
#include "G4LeptonConstructor.hh" 
#include "G4BaryonConstructor.hh" 
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


void ConstructParticle();

struct mapless 
{
	bool operator ()(const int i1,const int i2)
	{
	   return i1 < i2;
	}
};

void DoOneHadron(G4int maxEvents, G4int HadronPDG);

int main()
{
	G4int hadron;
	G4int maxEvents=1000;

// 	G4cout << " Welcome; How many events shal I do ? " ;
// 	G4cout.flush();
// 	G4cin >> maxEvents;
// 	
// 	G4cout << " Hadron PDG code " ;
// 	G4cout.flush();
// 	G4cin >> hadron;

	std::vector<int> hadrons;
	
	hadrons.push_back(2112); //proton
	hadrons.push_back(2112); //neutron
	hadrons.push_back(2224); // Delta++
	hadrons.push_back(2214); // Delta+
	hadrons.push_back(2114); // Delta0
	hadrons.push_back(1114); // Delta-
	
	hadrons.push_back(3122); // Lambda
	hadrons.push_back(3222); // Sigma+
	hadrons.push_back(3212); // Sigma0
	hadrons.push_back(3112); // Sigma-
	hadrons.push_back(3224); // Sigma*+
	hadrons.push_back(3214); // Sigma*0
	hadrons.push_back(3114); // Sigma*-
	hadrons.push_back(3322); // _-0
	hadrons.push_back(3312); // _--
	hadrons.push_back(3324); // _-*0
	hadrons.push_back(3314); // _-*-
	hadrons.push_back(3334); // Omega-
	
	

	ConstructParticle();
	
	std::vector<int>::iterator Ihadron;
	for (Ihadron=hadrons.begin();Ihadron!=hadrons.end();++Ihadron)
		DoOneHadron(maxEvents, *Ihadron);
	
}

	
void DoOneHadron(G4int maxEvents, G4int hadron)
{	
	std::map<const int, int, mapless> counts;				
//	std::pair<int,int> PDGcodes;
	const G4int PDGscale=10000;
	for (G4int i=0; i<maxEvents;  ++i)
	{
	    G4DiffractiveSplitableHadron diffHadron(
			G4ParticleTable::GetParticleTable()->FindParticle(hadron));
	    diffHadron.SplitUp();
	    G4int PDGcode1=diffHadron.GetNextParton()->GetPDGcode();
	    G4int PDGcode2=diffHadron.GetNextParton()->GetPDGcode();
	    int index=PDGscale*PDGcode1+PDGcode2;
	    counts[index]++;
// 	    G4cout << " test : count = " << counts[index] << G4endl;
// 	    G4cout << PDGcode1 << "  "<<
// 	    	      PDGcode2;
// 	    G4cout   << G4endl;	    
	}
	std::map<const int, int, mapless>::iterator iter;

	G4cout  << G4endl 
		<< G4ParticleTable::GetParticleTable()->
		                     FindParticle(hadron)->
				     GetParticleName() 
		<< " decays : " 
		<< G4endl;
	for (iter=counts.begin(); iter != counts.end();++iter) {
	
		G4cout << double(iter->second)/maxEvents*100 << "% "<<
//		iter->first/PDGscale << " - " << iter->first%PDGscale <<
			G4ParticleTable::GetParticleTable()->
		                     FindParticle(iter->first/PDGscale)->
				     GetParticleName() 
			<< " - " <<
			G4ParticleTable::GetParticleTable()->
		                     FindParticle(iter->first%PDGscale)->
				     GetParticleName()
			<< G4endl;
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
