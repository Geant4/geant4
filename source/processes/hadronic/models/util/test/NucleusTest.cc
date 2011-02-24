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
// $Id: NucleusTest.cc,v 1.4 2008-07-08 16:08:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 file
//
//
//             by Gunter Folger, June 1998.
//       class exercising G4Nucleus class.
//         testing static properties of all Isotopes
// ------------------------------------------------------------

#include "G4Fancy3DNucleus.hh"



#include "G4StableIsotopes.hh"
#include "G4NucleiPropertiesTable.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4ProtonField.hh"
#include "G4NeutronField.hh"

std::pair<G4double,G4double> AvergeMass(
		G4Fancy3DNucleus & nucleus, G4double A, G4double Z);

int main()
{

	G4Fancy3DNucleus nucleus;
	
// 	
// 	     
//
	
	G4StableIsotopes theIso;
		

	for (int Z=1; Z<93; Z++ )
	{
	   G4cout <<G4endl<< G4endl<< "new Element " ;
	   G4cout << theIso.GetName(Z) << G4endl;
	   G4cout << "         with " << theIso.GetNumberOfIsotopes(Z) ;
	   G4cout << " Isotopes"<< G4endl;
	   for (G4int iso=0; iso<theIso.GetNumberOfIsotopes(Z); iso++)
	   {
	       G4int massnumber=
	          theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+iso);
              for ( G4int repeat=0; repeat< 1; ++repeat )
	      {
	           nucleus.Init(massnumber,(G4double) Z);
		   nucleus.SortNucleonsInZ();
	      }
	      G4cout << "Charge" << nucleus.GetCharge() << G4endl;
	      G4cout << "Mass number" << nucleus.GetMassNumber() << G4endl;
	      G4cout << "nuclear (Table)Mass/Binding energy/ Ion mass: " << 
	        G4NucleiPropertiesTable::GetNuclearMass(Z,massnumber)
		<< " / " <<
		G4NucleiPropertiesTable::GetBindingEnergy(Z,massnumber)
		<< " / " <<
		G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,massnumber)
	        << G4endl;
	      std::pair<G4double,G4double> avmass=
	      			AvergeMass(nucleus, massnumber, (G4double) Z);
	      G4cout << "A, Z, N " << massnumber
	             << " " << Z
		     << " " << massnumber - Z
	      	     << " Mass " << nucleus.GetMass() 
	             << " Average mass " << avmass.first 
		     << " width " << avmass.second 
		     << " delta " << nucleus.GetMass()-avmass.first
		     << G4endl;
//	      G4cout << "Radius (.5) " << nucleus.GetNuclearRadius() << G4endl;
	      G4cout << "Radius (-2) " << nucleus.GetNuclearRadius(0.01) << G4endl;
//	      G4cout << "Radius (-4) " << nucleus.GetNuclearRadius(0.0001) << G4endl;
	      G4cout << "outermost nucleon " << nucleus.GetOuterRadius() << G4endl;

	      G4LorentzVector momentum;
	      G4Nucleon * anucleon;
	      nucleus.StartLoop();
	      while ( (anucleon=nucleus.GetNextNucleon() ) != NULL ) {
	             momentum = momentum + anucleon->GetMomentum();
	         }
	      G4cout << "Total Energy/Momentum (MeV) : " 
	      	     << momentum.e()/MeV << " " 
	      	     << momentum.x()/MeV << " "
		     << momentum.y()/MeV << " "
		     << momentum.z()/MeV << " "
		     << G4endl
		     << " Mass of sum(nucleons): " << momentum.mag() << " "
		     << G4endl;
	  
	  
	      G4bool print_positions=massnumber<10;
	        
	      if ( print_positions && nucleus.StartLoop() ) {
	         G4Nucleon * nucleon; 
	         G4cout << "radii" << G4endl;
	         while ( (nucleon=nucleus.GetNextNucleon() ) != NULL ) {
	            G4cout << "Mass/Energy/Momentum (MeV) : " 
	      	     << nucleon->GetMomentum().mag() << " " 
		     << nucleon->GetMomentum().e()/MeV << " " 
	      	     << nucleon->GetMomentum().x()/MeV << " "
		     << nucleon->GetMomentum().y()/MeV << " "
		     << nucleon->GetMomentum().z()/MeV << " "
		     << G4endl;
		 
	            G4cout << nucleon->GetPosition().mag() << G4endl;
	         }  
	      }
	   }   
	}
	return 0;
	
} 
std::pair<G4double,G4double> AvergeMass(
		G4Fancy3DNucleus & nucleus, G4double A, G4double Z)
{

	nucleus.Init(A,Z);
        G4ProtonField ProtonField(&nucleus);
	G4NeutronField NeutronField(&nucleus);
	G4int repeats=1000;
	G4double S_average=0, S_average2=0;
	for(G4int loop=0;loop<repeats; ++loop)
	{
	  nucleus.Init(A,Z);
	  G4LorentzVector momentum;
	  G4double fieldsum=0;	
	  G4Nucleon * anucleon;
	  nucleus.StartLoop();
	  while ( (anucleon=nucleus.GetNextNucleon() ) != NULL ) {
	  	  G4LorentzVector p_nuc=anucleon->GetMomentum();
		  G4double T=p_nuc.e() - p_nuc.mag();
		  p_nuc.setE(anucleon->GetDefinition()->GetPDGMass() + T);
		  p_nuc.setE(sqrt( sqr(anucleon->GetDefinition()->GetPDGMass()) 
		                    + p_nuc.vect().mag2()));
		  momentum += p_nuc;
//		  momentum += anucleon->GetMomentum();
		  G4double field=
		  (anucleon->GetParticleType() == G4Proton::Proton())
		  ? ProtonField.GetField(anucleon->GetPosition()) :
		    NeutronField.GetField(anucleon->GetPosition());
		  fieldsum += field;
	  }
	  S_average +=momentum.e() + fieldsum;            // field is negative
	  S_average2+=sqr(momentum.e() + fieldsum);
	// 		G4cout << "Total Energy(MeV) : " 
	// 	      		<< momentum.e()/MeV << " " 
	// 			<< G4endl;

	}
	std::pair<G4double,G4double> result;
	G4double average=S_average/repeats;
	result.first  =average;
	result.second =sqrt(std::fabs(
		repeats/(repeats-1) * (S_average2/repeats - sqr(average)) ));
	return result;
}
	
	
	
