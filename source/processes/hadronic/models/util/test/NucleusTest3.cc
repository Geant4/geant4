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
// $Id: NucleusTest3.cc,v 1.1 2003-10-08 13:48:24 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/HBookHistogram.h"
#include "CLHEP/Hist/HBookTuple.h"

void testNucleus(G4int maxNucleons, G4V3DNucleus & nucleus, G4int A, G4int Z);

#include "G4StableIsotopes.hh"
#include "G4FermiMomentum.hh"
#include "G4NuclearShellModelDensity.hh"
#include "G4NuclearFermiDensity.hh"


int main()
{
	HBookFile * hbfile=new HBookFile("Nucleustest2.hbook", 29); 
	G4Fancy3DNucleus nucleus;
	G4int maxNucleons;

	G4cout << " Welcome, how many Nucleons per Nucleus shall I do ? ";
	G4cin >> maxNucleons;
	G4cout << G4endl;

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

	hbfile->write();
	return 0;
}
	

void testNucleus(G4int maxNucleons, G4V3DNucleus & nucleus, G4int A, G4int Z)
{
	HBookTuple tuple(" theory density ",3000+A);
	
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
		tuple.column("radius",anR);
		tuple.column("density",
		          nucleardensity->GetDensity(G4ThreeVector(anR,0.,0.))*densityfactor);
		tuple.column("reldens",
		          nucleardensity->GetRelativeDensity(G4ThreeVector(anR,0.,0.)));
		tuple.dumpData();		          

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
	
	G4int momentumscale= maxFermi/MeV * 0.13;

	HBookHistogram density("Density distribution", 240,
	                            0., 12.*fermi, A); 
	HBookHistogram reldensity("relative Density distribution", 240,
	                            0., 12.*fermi, 500+A); 
	HBookHistogram momentum("Momentum distribution", 400,
	                            0., momentumscale * 10, 1000+A);
	                            
	HBookHistogram momentum_fraction(
	                      "Momentum distribution relativ",200,0., 1.2,
	                      2000+A);
	HBookHistogram momentum_fraction_last(
	                      "Momentum distribution relativ",200,0., 1.2,
	                      2500+A);
	HBookHistogram * momentum_frac[2];

	nucleus.Init((G4double) A,(G4double) Z);

	G4int charge=nucleus.GetCharge();
	G4cout << "Charge" << charge << G4endl;
	G4cout << "Mass number" << nucleus.GetMassNumber() << G4endl;
	G4cout << "Mass       " << nucleus.GetMass() << G4endl;
	G4cout << "Radius (-2) " << nucleus.GetNuclearRadius() << G4endl;

	G4int nNuclei = maxNucleons/A +1;

	for (G4int repeat=0; repeat< nNuclei; repeat++)
	{   
	    nucleus.Init((G4double) A,(G4double) Z);

//	       G4cout << "outer       " << nucleus.GetOuterRadius() << G4endl;

	    nucleus.StartLoop();
	    G4Nucleon * nucleon;
	    G4double radius,weight;
	    G4ThreeVector sum;
	    G4double MomentumFraction;
	    G4int nuc_count=0;
	    while ( (nucleon=nucleus.GetNextNucleon() ) != NULL ) {
		  radius = nucleon->GetPosition().mag();
		  weight=1./sqr(radius/fermi + 0.0001);
		  density.accumulate(radius);
		  reldensity.accumulate(radius,weight);
		  
		  G4double Pnucleon=nucleon->GetMomentum().vect().mag();
		 momentum.accumulate(Pnucleon);
		 sum+=Pnucleon;
		 if (Pnucleon< .0001 * MeV && A>1)
		 {
		    cout << G4endl << "iteration " << repeat << G4endl; 
		    G4cout << " got one very low momentum nucleon, event # "<< repeat << G4endl;
		    nucleus.StartLoop();
		    while ( (nucleon=nucleus.GetNextNucleon() ) != NULL ) {
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
		      cout << G4endl << "iteration " << repeat << G4endl; 
		      cout << " find error with fermi momentum" 
		           << " Nucleon/ here : (e, p, pfermi) "
		           << nucleon->Get4Momentum().e() << " / "
			   << nucleon->Get4Momentum().vect().mag() << " / "
		           << PnucFermi << G4endl;
		      cout << " Nucleon " << *nucleon << G4endl;
		 }
		 momentum_fraction.accumulate(MomentumFraction);
		     
	    }
	    momentum_fraction_last.accumulate(MomentumFraction);
	    if ( sum.mag() < perMillion*MeV && A>1) 
	    {
		G4cout << " event # "<< repeat << G4endl;
	    	G4cout << " Sum of momenta != 0 " << sum << sum.mag() << G4endl;
	    }
	}
	
	delete fermiMomentum;
	delete nucleardensity;
}
	
	
