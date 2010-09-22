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
// $Id: G4PreCompoundDeexcitation.cc,v 1.1 2010-09-22 22:17:08 yarba Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.

#include "G4PreCompoundDeexcitation.hh"

#include "G4PreCompoundModel.hh"

#include "globals.hh"

#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CascadeMomentum.hh"

// Constructor and destructor

G4PreCompoundDeexcitation::G4PreCompoundDeexcitation() 
  : G4CascadeColliderBase("G4PreCompoundDeexcitation"),
    theExcitationHandler( new G4ExcitationHandler )
{

   theDeExcitation = new G4PreCompoundModel( theExcitationHandler );

}

G4PreCompoundDeexcitation::~G4PreCompoundDeexcitation() 
{

   delete theExcitationHandler; // we need to delete here because G4PreComp does NOT delete it
   delete theDeExcitation;

}

// Main processing

void G4PreCompoundDeexcitation::collide(G4InuclParticle* /*bullet*/, 
			  	        G4InuclParticle* target,
				        G4CollisionOutput& globalOutput) 
{

  if (verboseLevel > 1) {
    G4cout << " >>> G4PreCompoundDeexcitation::collide" << G4endl;
  }
  
  // Ensure that input state is sensible
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(target);
  if (!ntarget) {
    G4cerr << " G4PreCompoundDeexcitation ERROR:  residual fragment must be G4InuclNuclei"
	   << G4endl;
    return;
  }
  
  std::vector<G4InuclElementaryParticle> particles;
  
  if ( ntarget->getExitationEnergy() > 0. && ntarget->getA() > 1.5 )
  {
      getDeExcitedFragments( ntarget, particles );  
  }
  else
  {
     convertFragment( ntarget, particles );
  }
   
  globalOutput.addOutgoingParticles(particles);		// Evaporated particles and nucleus
  
  return;
  
}
  
void G4PreCompoundDeexcitation::convertFragment( G4InuclNuclei* rfrag, 
                                                 std::vector<G4InuclElementaryParticle> particles )
{

   G4ParticleTable* theTableOfParticles = G4ParticleTable::GetParticleTable();

   const G4CascadeMomentum& mom = rfrag->getMomentum();
   G4int A = G4int(rfrag->getA());
   G4int Z = G4int(rfrag->getZ());
   G4double eKin = rfrag->getKineticEnergy() * GeV;
   G4ParticleDefinition * aIonDef = theTableOfParticles->FindIon(Z, A, 0, Z);
   G4ThreeVector aMom(mom[1], mom[2], mom[3]);
   aMom = aMom.unit();
   G4DynamicParticle aFragment( aIonDef, aMom, eKin );
   particles.push_back( G4InuclElementaryParticle(aFragment) );

   return;
   
}

void G4PreCompoundDeexcitation::getDeExcitedFragments( G4InuclNuclei* rfrag, 
                                                       std::vector<G4InuclElementaryParticle> particles )
{

  const G4CascadeMomentum& mom = rfrag->getMomentum();
  G4int A = G4int(rfrag->getA());
  G4int Z = G4int(rfrag->getZ());
  G4LorentzVector aMomentum = G4LorentzVector(mom[1]*GeV,mom[2]*GeV,
					      mom[3]*GeV,mom[0]*GeV);
  G4ExitonConfiguration exiton = rfrag->getExitonConfiguration();
    
  G4Fragment frag(A,Z,aMomentum);
  frag.SetParticleDefinition( (rfrag->getDynamicParticle()).GetDefinition() );
  
  // this is a dummy method, the exec complains about it at run time
  //
  // frag.SetExcitationEnergy((rfrag->getExitationEnergy())*MeV);
  
  
  frag.SetNumberOfHoles((G4int)(exiton.protonHoles+exiton.neutronHoles));
  frag.SetNumberOfParticles((G4int)(exiton.protonQuasiParticles+exiton.protonQuasiParticles));
  frag.SetNumberOfCharged((G4int)(exiton.protonQuasiParticles));

  G4ReactionProductVector* precompoundProducts = 0;
  
  if (  ( explosion(rfrag) || Z==0 ) && theExcitationHandler ) // in principle, the explosion(...) stuff should also 
                                                     // handle properly the case of Z=0 (neutron blob)
  {
     precompoundProducts=theExcitationHandler->BreakItUp(frag);
  }
  else
  {
     precompoundProducts = theDeExcitation->DeExcite(frag);
  }

  //  G4LorentzVector pSumPreco(0), pPreco(0);
  if ( precompoundProducts )  {
    std::vector<G4ReactionProduct *>::iterator j;
    for(j=precompoundProducts->begin(); j!=precompoundProducts->end(); ++j) {
      G4DynamicParticle aFragment((*j)->GetDefinition(), 
				  (*j)->GetMomentum(), 
				  (*j)->GetKineticEnergy());

      particles.push_back( G4InuclElementaryParticle(aFragment) );
    }
    precompoundProducts->clear();
    delete precompoundProducts;
  }

   return;
}
