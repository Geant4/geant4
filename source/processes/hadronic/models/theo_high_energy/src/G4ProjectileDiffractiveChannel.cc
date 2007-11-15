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
// $Id: G4ProjectileDiffractiveChannel.cc,v 1.2 2007-11-15 16:07:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Author : Gunter Folger Nov 2007
// Class Description
// Final state production model for theoretical models of hadron inelastic
// projectile diffractive scattering in geant4;
// Class Description - End
//
// Modified:

#include "G4ProjectileDiffractiveChannel.hh"

#include "G4HadTmpUtil.hh"		//lrint
#include "G4QHadronVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Lambda.hh"

//#define debug_getFraction 1
//#define debug_scatter 1

G4ProjectileDiffractiveChannel::G4ProjectileDiffractiveChannel()
{
	theQDiffraction=G4QDiffractionRatio::GetPointer();
}

G4ProjectileDiffractiveChannel::~G4ProjectileDiffractiveChannel()
{}

G4double G4ProjectileDiffractiveChannel::GetFraction(G4Nucleus &theNucleus,
				const G4DynamicParticle & thePrimary)
{
	G4double ratio;
	ratio=theQDiffraction->GetRatio(thePrimary.GetTotalMomentum(),
			thePrimary.GetDefinition()->GetPDGEncoding(),
			G4lrint(theNucleus.GetZ()),
			G4lrint(theNucleus.GetN()-theNucleus.GetZ()));
	   #ifdef debug_getFraction			
	      G4cout << "G4ProjectileDiffractiveChannel::ratio " << ratio << G4endl;
	   #endif
	       
	return ratio;       		
}

G4KineticTrackVector * G4ProjectileDiffractiveChannel::Scatter(G4Nucleus &theNucleus,
					const G4DynamicParticle & thePrimary)
{
	
	
	G4int A=G4lrint(theNucleus.GetN());
	G4int Z=G4lrint(theNucleus.GetZ());

	G4ParticleDefinition * pDef=thePrimary.GetDefinition();
	   #ifdef debug_scatter
	     G4cout << " ProjectileDiffractive: A, Z,  proj-pdg" <<" "<< A <<" "<<Z 
			<< " "<<" " << pDef->GetParticleName()<< G4endl;
	   #endif
	
	G4QHadronVector * result(0);
	G4KineticTrackVector * ktv(0);
	G4bool tryAgain;
	
	do {
	   tryAgain = false;
	   result=theQDiffraction->ProjFragment(pDef->GetPDGEncoding(),
			        thePrimary.Get4Momentum(),Z, A-Z);
				
	   if ( result && result->size() > 0 )
	   {
	       ktv=new G4KineticTrackVector();
	       std::vector<G4QHadron *>::iterator i;
		 #ifdef debug_scatter
		    G4int count(0);
		 #endif
	       for (i=result->begin(); i!=result->end(); i++)
	       {
		  G4int pdgCode=(*i)->GetPDGCode();
		  G4ParticleDefinition * secDef(0);
		  if ( pdgCode < 80000000) {     // check if ion; should be 90Million, but 89Million is possible
		      secDef= G4ParticleTable::GetParticleTable()->FindParticle(pdgCode);
		  }else {
		      G4int N = pdgCode % 1000;
		      G4int Z = (pdgCode/1000) %1000;
		      if ( Z < 500 && N < 500){   // protect for delta being coded as/in nucleus 
			 secDef = G4ParticleTable::GetParticleTable()->GetIon(Z,N+Z, 0);
			 if ( ! secDef ) 
			 {  // exceptions to the rule!
	        	    if ( pdgCode == 90000001 ) secDef= G4Neutron::Neutron();
			    if ( pdgCode == 91000000 ) secDef= G4Lambda::Lambda();
			 }
		      }	  
   		  }
		  
        	  if  ( secDef ){

		      G4ThreeVector pos=(*i)->GetPosition();  //GetPosition returns const &
		      G4LorentzVector mom=(*i)->Get4Momentum();

		      G4KineticTrack * sPrim=new G4KineticTrack(secDef,
				      (*i)->GetFormationTime(), pos, mom);
		      ktv->push_back(sPrim);

			#ifdef debug_scatter
			  G4cout << "G4ProjectileDiffractive sec # "  << ++count << ", "
	        	    << "ChipsPDGCode=" << pdgCode << ", "
	        	    << secDef->GetParticleName() 
	        	    << ", time(ns) " << (*i)->GetFormationTime()/ns
	        	    << ", pos " << pos
	        	    << ", 4mom " << (*i)->Get4Momentum()
	        	    << G4endl;
			#endif
	           } else {
		       #ifdef debug_scatter
		         G4cout << "G4ProjectileDiffractiveChannel: Particle with PDG code "<< pdgCode <<" does not converted!!!"<<G4endl;
		       #endif 
		       tryAgain=true;
		   }
	       }   
	       std::for_each(result->begin(), result->end(), DeleteQHadron());
               delete result;
	       result=0;
	   }
	    
	   if ( result ) delete result;
	   if ( tryAgain ) {     // QDiffraction returned a "difficult" particle in nucleus
	       std::for_each(ktv->begin(), ktv->end(), DeleteKineticTrack());
	       delete ktv;
	       ktv=0;
	   }
	} while ( ! ktv);          

	return ktv;			       
}
