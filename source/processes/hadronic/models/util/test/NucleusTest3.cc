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
// $Id$
//
// ------------------------------------------------------------
//      GEANT 4 file
//
//
//             by Gunter Folger, June 1998.
//       class exercising G4Nucleus class.
//         look at density and momentum distributions.
// ------------------------------------------------------------



//-------------------------------------------------------------
#include "G4Fancy3DNucleus.hh"


#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

void testNucleus(G4int maxNucleons, G4V3DNucleus & nucleus, G4int A, G4int Z);

#include "G4StableIsotopes.hh"
#include "G4FermiMomentum.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NuclearFermiDensity.hh"

#include <iostream>
#include <sstream>


int main()
{
	G4Fancy3DNucleus nucleus;
	G4int maxNucleons;

	G4cout << " Welcome, how many Nucleons per Nucleus shall I do ? ";
	G4cin >> maxNucleons;
	G4cout << G4endl;
	TFile * rootFile = new TFile("NucleusTest3.root","CREATE");
	if ( ! rootFile ) 
	{
	   std::cout << " Fail to create root file " << std::endl;
	   exit(1);
	}   

//	testNucleus(maxNucleons, nucleus,1,1);
	testNucleus(maxNucleons, nucleus,2,1);
	testNucleus(maxNucleons, nucleus,3,2);
	testNucleus(maxNucleons, nucleus,4,2);
	
	
	testNucleus(maxNucleons, nucleus,16,8);
	
	testNucleus(maxNucleons, nucleus,17,8);
	
	testNucleus(maxNucleons, nucleus,56,26);
	
	testNucleus(maxNucleons, nucleus,137,55);
	 
	testNucleus(maxNucleons, nucleus,197,79); 
	 
	testNucleus(maxNucleons, nucleus,207,82); 

	rootFile->Write();
	return 0;
}
	

void testNucleus(G4int maxNucleons, G4V3DNucleus & nucleus, G4int A, G4int Z)
{
        static bool init(true);
	static std::vector<TTree *> tree_vector;
	static std::vector<TH1F *>  hist_vector;
	static TTree * Tdensity;
        struct theo_density {
	     int A;
	     float radius;	     
	     float rho;
	     float relrho;
	};
	static theo_density t_density;     
	if ( init) {
	   init=false;

   //	std::cout << " str-A : " << str_A.str() << std::endl;
	   TString TTitle("Rho");
           Tdensity=new TTree("Rho","density");
	   Tdensity->Branch("theodens",&t_density.A,"massnum/I:radius/F:rho/F:relrho/F");

	   tree_vector.push_back(Tdensity);
	
	}
	G4StableIsotopes theIso;
	G4VNuclearDensity * nucleardensity;
	G4FermiMomentum * fermiMomentum;
	
	if ( A < 17 ) {
	    nucleardensity = new G4NuclearShellModelDensity(A, Z);
	} else {
	    nucleardensity = new G4NuclearFermiDensity(A, Z);
	}
	
	G4double densityfactor= A * pow(fermi,3);
	for (G4double anR=0.; anR< 10.*fermi; anR+= 0.001*fermi)
	{
		t_density.A=A;
		t_density.radius=anR;
		t_density.rho=nucleardensity->GetDensity(G4ThreeVector(anR,0.,0.))*densityfactor;
		t_density.relrho=nucleardensity->GetRelativeDensity(G4ThreeVector(anR,0.,0.));
		Tdensity->Fill();

	}
	fermiMomentum= new G4FermiMomentum();
	fermiMomentum->Init(A,Z);
	
	G4double radiusforMaxDensity=nucleardensity->GetRadius(1.);
        
	G4double maxDensity=nucleardensity->GetDensity(G4ThreeVector(radiusforMaxDensity, 0.,0.));

	G4double maxFermi = fermiMomentum->GetFermiMomentum(maxDensity);
	
	G4cout <<G4endl<< G4endl<< "new Element " ;
	G4cout << theIso.GetName(Z) << G4endl;
	
	G4cout << "relative density = 1 at r= " << radiusforMaxDensity << G4endl;
	G4cout << " Max density: " << A*maxDensity * pow(fermi,3)
		<< "( Nucleons/fermi^3 )" << G4endl;
	G4cout << " density(R=0): " 
		<< A*nucleardensity->GetDensity(G4ThreeVector(0., 0., 0.))* pow(fermi,3)
		<< "( Nucleons/fermi^3 )" << G4endl;

	G4cout << " Max Fermi Momentum " << maxFermi /MeV<< " MeV"<< G4endl;;
	
	G4int momentumscale= G4int(maxFermi/MeV * 0.13);

        std::stringstream str_A(std::stringstream::out);
        str_A << "A=" << A;
	TString title("density for "+str_A.str());
	std::cout << " root title: " << title << " str_A: " << str_A.str() << std::endl;
	const char * T_A = str_A.str().data();
	TH1F * density=new TH1F(("Density distribution "+str_A.str()).data(),T_A, 
	     		     240, 0., 12.*fermi);
	hist_vector.push_back(density); 
	TH1F * reldensity=new TH1F(("relative Density distribution "+str_A.str()).data(),T_A, 
	     		     240, 0., 12.*fermi ); 
	hist_vector.push_back(reldensity); 
	TH1F * momentum=new TH1F(("Momentum distribution "+str_A.str()).data(),T_A, 
	     		     400, 0., momentumscale * 10);
	     			   
	hist_vector.push_back(momentum); 
	TH1F * momentum_fraction=new TH1F(("Momentum distribution relativ "+str_A.str()).data(),T_A,
	        	     200,0., 1.2);
	hist_vector.push_back(momentum_fraction); 
	TH1F * momentum_fraction_last=new TH1F(("Momentum distribution relativ "+str_A.str()).data(),T_A,
	                      200,0., 1.2);
	hist_vector.push_back(momentum_fraction_last); 
//	THF1 * momentum_frac[2];

	nucleus.Init( A, Z);

	G4int charge=nucleus.GetCharge();
	G4cout << "Charge" << charge << G4endl;
	G4cout << "Mass number" << nucleus.GetMassNumber() << G4endl;
	G4cout << "Mass       " << nucleus.GetMass() << G4endl;
	G4cout << "Radius (-2) " << nucleus.GetNuclearRadius() << G4endl;

	G4int nNuclei = maxNucleons/A +1;
        G4double norm=1./nNuclei;
	for (G4int repeat=0; repeat< nNuclei; repeat++)
	{   
	    nucleus.Init( A, Z);

//	       G4cout << "outer       " << nucleus.GetOuterRadius() << G4endl;

	    nucleus.StartLoop();
	    G4Nucleon * nucleon;
	    G4double radius,weight;
	    G4ThreeVector sum;
	    G4double MomentumFraction;
	    G4int nuc_count=0;
	    while ( (nucleon=nucleus.GetNextNucleon()) ) {
		  radius = nucleon->GetPosition().mag();
		  weight=1./sqr(radius/fermi + 0.0001);
		  density->Fill(radius,norm);
		  reldensity->Fill(radius,norm*weight);
		  
		  G4double Pnucleon=nucleon->GetMomentum().vect().mag();
		 momentum->Fill(Pnucleon);
		 sum+=nucleon->GetMomentum().vect();
		 if (Pnucleon< .0001 * MeV && A>1)
		 {
		    G4cout << G4endl << "iteration " << repeat << G4endl; 
		    G4cout << " got one very low momentum nucleon, event # "<< repeat << G4endl;
		    nucleus.StartLoop();
		    while ( (nucleon=nucleus.GetNextNucleon() ) ) {
		    G4cout << " Position " << nucleon->GetPosition() << G4endl;
		    G4cout << " Momentum " << nucleon->GetMomentum() << G4endl;
		    }
		   break; 
		 }
		 G4double PnucFermi=
		       fermiMomentum->GetFermiMomentum(nucleardensity->
		 	   GetDensity(nucleon->GetPosition()));
		 MomentumFraction=Pnucleon/PnucFermi;
		 if ( PnucFermi < nucleon->Get4Momentum().vect().mag() )
		 {
		      G4cout << G4endl << "iteration " << repeat << G4endl; 
		      G4cout << " find error with fermi momentum" 
		           << " Nucleon/ here : (e, p, pfermi) "
		           << nucleon->Get4Momentum().e() << " / "
			   << nucleon->Get4Momentum().vect().mag() << " / "
		           << PnucFermi << G4endl;
		      G4cout << " Nucleon " << *nucleon << G4endl;
		 }
		 momentum_fraction->Fill(MomentumFraction);
		     
	    }
	    momentum_fraction_last->Fill(MomentumFraction);
	    if ( sum.mag() > perMillion*MeV && A>1) 
	    {
		G4cout << " event # "<< repeat << G4endl;
	    	G4cout << " Sum of momenta != 0 " << sum << sum.mag() << G4endl;
	    }
	}
	
	delete fermiMomentum;
	delete nucleardensity;
}
	
	
