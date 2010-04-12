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
// $Id: G4InuclEvaporation.cc,v 1.13 2010-04-12 23:39:41 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100405  M. Kelsey -- Pass const-ref std::vector<>
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide(), use
//		const_iterator.

#include <numeric>
#include "G4IonTable.hh"
#include "globals.hh"
#include "G4V3DNucleus.hh"
#include "G4DynamicParticleVector.hh"
#include "G4EvaporationInuclCollider.hh"
#include "G4InuclEvaporation.hh"
#include "G4InuclNuclei.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4Nucleon.hh"
#include "G4NucleiModel.hh"
#include "G4HadronicException.hh"
#include "G4LorentzVector.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"

typedef std::vector<G4InuclElementaryParticle>::const_iterator particleIterator;
typedef std::vector<G4InuclNuclei>::const_iterator nucleiIterator;


G4InuclEvaporation::G4InuclEvaporation() {
  verboseLevel=0;
}

G4InuclEvaporation::G4InuclEvaporation(const G4InuclEvaporation &) : G4VEvaporation() {
  throw G4HadronicException(__FILE__, __LINE__, "G4InuclEvaporation::copy_constructor meant to not be accessable.");
}


G4InuclEvaporation::~G4InuclEvaporation() {
}

const G4InuclEvaporation & G4InuclEvaporation::operator=(const G4InuclEvaporation &) {
  throw G4HadronicException(__FILE__, __LINE__, "G4InuclEvaporation::operator= meant to not be accessable.");
  return *this;
}


G4bool G4InuclEvaporation::operator==(const G4InuclEvaporation &) const {
  return false;
}

G4bool G4InuclEvaporation::operator!=(const G4InuclEvaporation &) const {
  return true;
}

void G4InuclEvaporation::setVerboseLevel( const G4int verbose ) {
  verboseLevel = verbose;
}

G4FragmentVector * G4InuclEvaporation::BreakItUp(const G4Fragment &theNucleus) {
 
  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3,
                      pionMinus = 5, pionZero = 7, photon = 10,
                      kaonPlus = 11, kaonMinus = 13, kaonZero = 15,
                      kaonZeroBar = 17, lambda = 21, sigmaPlus = 23,
                      sigmaZero = 25, sigmaMinus = 27, xiZero = 29, xiMinus = 31 };  

  std::vector< G4DynamicParticle * > secondaryParticleVector;
  G4FragmentVector * theResult = new G4FragmentVector;

  if (theNucleus.GetExcitationEnergy() <= 0.0) { // Check that Excitation Energy > 0
    theResult->push_back(new G4Fragment(theNucleus));
    return theResult;
  }

  G4double A = theNucleus.GetA();
  G4double Z = theNucleus.GetZ();
  G4double mTar  = G4NucleiProperties::GetNuclearMass(A, Z); // Mass of the target nucleus
  G4LorentzVector tmp =theNucleus.GetMomentum();

  G4ThreeVector momentum = tmp.vect();
  //  G4double energy = tmp.e();
  G4double exitationE = theNucleus.GetExcitationEnergy();

  // Move to CMS frame, save initial velocity of the nucleus to boostToLab vector.
  //   G4ThreeVector boostToLab( ( 1/G4NucleiProperties::GetNuclearMass( A, Z ) ) * momentum ); 
  G4InuclNuclei* tempNuc = new G4InuclNuclei(A, Z);
  G4double mass=tempNuc->getMass()*1000;
  G4ThreeVector boostToLab( ( 1/mass) * momentum ); 

  if ( verboseLevel > 2 )
    G4cout << " G4InuclEvaporation : initial kinematics : boostToLab vector = " << boostToLab << G4endl
	   << "                     excitation energy  : " << exitationE << G4endl;

  if (verboseLevel > 2) {
    G4cout << "G4InuclEvaporation::BreakItUp >>> A: " << A << " Z: " << Z << " exitation E: " << 
      exitationE << " mass: " << mTar << G4endl;
  };

  G4InuclNuclei* nucleus = new G4InuclNuclei(A, Z);
  nucleus->setExitationEnergy(exitationE/1000);

  G4EquilibriumEvaporator*          eqil = new G4EquilibriumEvaporator;
  G4Fissioner*                      fiss = new G4Fissioner;
  G4BigBanger*                      bigb = new G4BigBanger;
  G4EvaporationInuclCollider* evaporator = new G4EvaporationInuclCollider(eqil, fiss, bigb);

  G4CollisionOutput output;
  evaporator->collide(0, nucleus, output);

  const std::vector<G4InuclNuclei>& nucleiFragments = output.getNucleiFragments();
  const std::vector<G4InuclElementaryParticle>& particles = output.getOutgoingParticles();

  G4double ekin,emas;
  G4double eTot=0.0;
  G4DynamicParticle* cascadeParticle = 0;
  G4int  i=1;
  //G4cout << "# particles: " << output.getOutgoingParticles().size() << G4endl;
  if (!particles.empty()) { 
    particleIterator ipart;
    G4int outgoingParticle;

    for (ipart = particles.begin(); ipart != particles.end(); ipart++) {
     
      outgoingParticle = ipart->type();

      if (verboseLevel > 2) {
	G4cout << "Evaporated particle:  " << i << " of type: " << outgoingParticle << G4endl;
        i++;
	//       	ipart->printParticle();
      }

      eTot += ipart->getEnergy()*1000;
      ekin = ipart->getKineticEnergy()*1000;
      emas = ipart->getMass()*1000;

      G4LorentzVector v(ipart->getMomentum().vect()*1000, (ekin+emas));
      v.boost( boostToLab );

      switch(outgoingParticle) {

      case proton: 
	cascadeParticle = new G4DynamicParticle(G4Proton::ProtonDefinition(), v.vect(), v.e());
	break; 

      case neutron: 
	cascadeParticle = new G4DynamicParticle(G4Neutron::NeutronDefinition(), v.vect(), v.e());
	break;

      case photon: 
	cascadeParticle = new G4DynamicParticle(G4Gamma::Gamma(), v.vect(), v.e());
	break;
      default:
        G4cout << " ERROR: GInuclEvapration::Propagate undefined particle type" << G4endl;
      };

      secondaryParticleVector.push_back( cascadeParticle );
      //      theResult.AddSecondary(cascadeParticle); 
    }  
  }
  fillResult( secondaryParticleVector, theResult);
  

  //  G4cout << "# fragments " << output.getNucleiFragments().size() << G4endl;
  i=1; 
  if (!nucleiFragments.empty()) { 
    nucleiIterator ifrag;

    for (ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) {

      ekin = ifrag->getKineticEnergy()*1000;
      emas = ifrag->getMass()*1000;
      eTot += ifrag->getEnergy()*1000;

      G4LorentzVector v(ifrag->getMomentum().vect()*1000, (ekin+emas));
      v.boost( boostToLab );

      if (verboseLevel > 2) {
	G4cout << " Nuclei fragment: " << i << G4endl; i++;
      }

      G4int A = G4int(ifrag->getA());
      G4int Z = G4int(ifrag->getZ());
      if (verboseLevel > 2) {
	G4cout << "boosted v" << v << G4endl;
      }
      theResult->push_back( new G4Fragment(A, Z, v) ); 
    }
  }

  //G4cout << ">>>> G4InuclEvaporation::BreakItUp end " << G4endl;
  return theResult;
}

void G4InuclEvaporation::fillResult(const std::vector<G4DynamicParticle*>& secondaryParticleVector,
				     G4FragmentVector * aResult )
{
  // Fill the vector pParticleChange with secondary particles stored in vector.
  for ( size_t i = 0 ; i < secondaryParticleVector.size() ; i++ )
    {
      G4int aZ = static_cast<G4int> (secondaryParticleVector[i]->GetDefinition()->GetPDGCharge() );
      G4int aA = static_cast<G4int> (secondaryParticleVector[i]->GetDefinition()->GetBaryonNumber());
      G4LorentzVector aMomentum = secondaryParticleVector[i]->Get4Momentum();
      if(aA>0) {
	aResult->push_back( new G4Fragment(aA, aZ, aMomentum) ); 
      } else {
	aResult->push_back( new G4Fragment(aMomentum, secondaryParticleVector[i]->GetDefinition()) ); 
      }
    }
  return;
}



