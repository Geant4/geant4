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
// $Id: G4AblaEvaporation.cc,v 1.1 2008-02-27 18:31:11 miheikki Exp $
//
#include <numeric>
// #include "G4IonTable.hh"
// #include "globals.hh"
// #include "G4V3DNucleus.hh"
// #include "G4DynamicParticleVector.hh"
// #include "G4EvaporationInuclCollider.hh"
// #include "G4InuclEvaporation.hh"
// #include "G4InuclNuclei.hh"
// #include "G4Track.hh"
// #include "G4Nucleus.hh"
// #include "G4Nucleon.hh"
// #include "G4NucleiModel.hh"
#include "G4HadronicException.hh"
// #include "G4LorentzVector.hh"
// #include "G4EquilibriumEvaporator.hh"
// #include "G4Fissioner.hh"
// #include "G4BigBanger.hh"
// #include "G4InuclElementaryParticle.hh"
// #include "G4InuclParticle.hh"
// #include "G4CollisionOutput.hh"

#include "G4AblaEvaporation.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"

G4AblaEvaporation::G4AblaEvaporation() {
  verboseLevel=0;
 hazard = new G4Hazard();
  // set initial values:
  // First random seed:
  // (Premiere graine)
  //  hazard->ial = 38035;
  hazard->ial = 979678188;
  // other seeds:
  hazard->igraine[0] = 3997;
  hazard->igraine[1] = 15573;
  hazard->igraine[2] = 9971;
  hazard->igraine[3] = 9821; 
  hazard->igraine[4] = 99233; 
  hazard->igraine[5] = 11167; 
  hazard->igraine[6] = 12399;
  hazard->igraine[7] = 11321; 
  hazard->igraine[8] = 9825;
  hazard->igraine[9] = 2587; 
  hazard->igraine[10] = 1775;
  hazard->igraine[11] = 56799; 
  hazard->igraine[12] = 1156;
  //  hazard->igraine[13] = 11207;
  hazard->igraine[13] = 38957; 
  hazard->igraine[14] = 35779; 
  hazard->igraine[15] = 10055; 
  hazard->igraine[16] = 76533; 
  hazard->igraine[17] = 33759;
  hazard->igraine[18] = 13227;
}

G4AblaEvaporation::G4AblaEvaporation(const G4AblaEvaporation &) : G4VEvaporation() {
  throw G4HadronicException(__FILE__, __LINE__, "G4AblaEvaporation::copy_constructor meant to not be accessable.");
}

G4AblaEvaporation::~G4AblaEvaporation() {
}

const G4AblaEvaporation & G4AblaEvaporation::operator=(const G4AblaEvaporation &) {
  throw G4HadronicException(__FILE__, __LINE__, "G4AblaEvaporation::operator= meant to not be accessable.");
  return *this;
}

G4bool G4AblaEvaporation::operator==(const G4AblaEvaporation &) const {
  return false;
}

G4bool G4AblaEvaporation::operator!=(const G4AblaEvaporation &) const {
  return true;
}

void G4AblaEvaporation::setVerboseLevel( const G4int verbose ) {
  verboseLevel = verbose;
}

G4FragmentVector * G4AblaEvaporation::BreakItUp(const G4Fragment &theNucleus) {
 

  G4VarNtp *varntp = new G4VarNtp();
  G4Volant *volant = new G4Volant();

  G4Abla *abla = new G4Abla(hazard, volant, varntp);
  G4cout <<"Initializing evaporation..." << G4endl;
  abla->initEvapora();
  G4cout <<"Initialization complete!" << G4endl;
  
  G4double nucleusA = theNucleus.GetA();
  G4double nucleusZ = theNucleus.GetZ();
  G4double nucleusMass = G4NucleiProperties::GetAtomicMass(nucleusA, nucleusZ);
  G4double excitationEnergy = theNucleus.GetExcitationEnergy();
  G4double angularMomentum = 0.0; // Don't know how to get this quantity... From Geant4???

  G4LorentzVector tmp = theNucleus.GetMomentum();

  G4ThreeVector momentum = tmp.vect();

  G4double recoilEnergy = tmp.e();
  G4double momX = momentum.x();
  G4double momY = momentum.y();
  G4double momZ = momentum.z();
  //  G4double energy = tmp.e();
  G4double exitationE = theNucleus.GetExcitationEnergy() * MeV;

  varntp->ntrack = -1;
  varntp->massini = theNucleus.GetA();
  varntp->mzini = theNucleus.GetZ();

  std::vector<G4DynamicParticle*> cascadeParticles;
  G4FragmentVector * theResult = new G4FragmentVector;
  if (theNucleus.GetExcitationEnergy() <= 0.0) { // Check that Excitation Energy > 0
    theResult->push_back(new G4Fragment(theNucleus));
    return theResult;
  }

  //  G4double mTar  = G4NucleiProperties::GetAtomicMass(A, Z); // Mass of the target nucleus
  varntp->exini = exitationE;

  G4int particleI, n = 0;

  // Print diagnostic messages. 0 = silent, 1 and 2 = verbose
  //  verboseLevel = 3;

  // Increase the event number:
  eventNumber++;

  G4DynamicParticle *cascadeParticle = 0;
  //  G4ParticleDefinition *aParticleDefinition = 0;

  // Map Geant4 particle types to corresponding INCL4 types.
  enum bulletParticleType {nucleus = 0, proton = 1, neutron = 2, pionPlus = 3, pionZero = 4, 
                           pionMinus = 5, deuteron = 6, triton = 7, he3 = 8, he4 = 9};

  // Check wheter the input is acceptable. This will contain more tests in the future. 

//   void breakItUp(G4double nucleusA, G4double nucleusZ, G4double nucleusMass, G4double excitationEnergy,
// 		       G4double angularMomentum, G4double recoilEnergy, G4double momX, G4double momY, G4double momZ)
  G4cout <<"Calling the actual ABLA model..." << G4endl;
  G4cout <<"Excitation energy: " << excitationEnergy << G4endl;
  abla->breakItUp(nucleusA, nucleusZ, nucleusMass, excitationEnergy, angularMomentum, recoilEnergy, momX, momY, momZ,
		  eventNumber);
  G4cout <<"Done." << G4endl;

  if(verboseLevel > 0) {
    // Diagnostic output
    G4cout <<"G4AblaEvaporation: Target A:  " << nucleusA << G4endl;
    G4cout <<"G4AblaEvaporation: Target Z:  " << nucleusZ << G4endl;

    for(particleI = 0; particleI < varntp->ntrack; particleI++) {
      G4cout << n << " ";
      G4cout << varntp->massini << " " << varntp->mzini << " ";
      G4cout << varntp->exini << " " << varntp->mulncasc << " " << varntp->mulnevap << " " << varntp->mulntot << " ";
      G4cout << varntp->bimpact << " " << varntp->jremn << " " << varntp->kfis << " " << varntp->estfis << " ";
      G4cout << varntp->izfis << " " << varntp->iafis << " " << varntp->ntrack << " " << varntp->itypcasc[particleI] << " ";
      G4cout << varntp->avv[particleI] << " " << varntp->zvv[particleI] << " " << varntp->enerj[particleI] << " ";
      G4cout << varntp->plab[particleI] << " " << varntp->tetlab[particleI] << " " << varntp->philab[particleI] << G4endl;
    }
  }

  // Loop through the INCL4+ABLA output.
  G4double momx, momy, momz; // Momentum components of the outcoming particles.
  G4double eKin;
  G4cout <<"varntp->ntrack = " << varntp->ntrack << G4endl;
  for(particleI = 0; particleI < varntp->ntrack; particleI++) {
    // Get energy/momentum and construct momentum vector:
    // In INCL4 coordinates!
    momx = varntp->plab[particleI]*std::cos(varntp->tetlab[particleI]*CLHEP::pi/180.0)*std::sin(varntp->philab[particleI]*CLHEP::pi/180.0)*MeV;
    momy = varntp->plab[particleI]*std::sin(varntp->tetlab[particleI]*CLHEP::pi/180.0)*std::sin(varntp->philab[particleI]*CLHEP::pi/180.0)*MeV;
    momz = varntp->plab[particleI]*std::cos(varntp->tetlab[particleI]*CLHEP::pi/180.0)*MeV;

    eKin = varntp->enerj[particleI] * MeV;

    if(verboseLevel > 1) {
      //      G4cout <<"Momentum direction: (x ,y,z)";
      //      G4cout << "(" << momx <<"," << momy << "," << momz << ")" << G4endl;
    }

    // This vector tells the direction of the particle.
    G4ThreeVector momDirection(momx, momy, momz);
    momDirection = momDirection.unit();
        
    // Identify the particle/nucleus:
    G4int particleIdentified = 0;

    // Proton
    if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 1)) {
      cascadeParticle = 
	new G4DynamicParticle(G4Proton::ProtonDefinition(), momDirection, eKin);
      particleIdentified++;
    }

    // Neutron
    if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 0)) {
      cascadeParticle = 
	new G4DynamicParticle(G4Neutron::NeutronDefinition(), momDirection, eKin);
      particleIdentified++;
    }

    // PionPlus
    if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 1)) {
      cascadeParticle = 
	new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), momDirection, eKin);
      particleIdentified++;
    }

    // PionZero
    if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 0)) {
      cascadeParticle = 
	new G4DynamicParticle(G4PionZero::PionZeroDefinition(), momDirection, eKin);
      particleIdentified++;
    }

    // PionMinus
    if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == -1)) {
      cascadeParticle = 
	new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), momDirection, eKin);
      particleIdentified++;
    }

    // Nuclei fragment
    if((varntp->avv[particleI] > 1) && (varntp->zvv[particleI] >= 1)) {
      G4ParticleDefinition * aIonDef = 0;
      G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();

      G4int A = G4int(varntp->avv[particleI]);
      G4int Z = G4int(varntp->zvv[particleI]);
      aIonDef = theTableOfParticles->FindIon(Z, A, 0, Z);
	
      cascadeParticle = 
	new G4DynamicParticle(aIonDef, momDirection, eKin);
      particleIdentified++;
    }

    // Check that the particle was identified properly.
    if(particleIdentified == 1) {
      // Put data into G4HadFinalState:
      cascadeParticle->Set4Momentum(cascadeParticle->Get4Momentum());
      cascadeParticles.push_back(cascadeParticle);
      //      theResult.AddSecondary(cascadeParticle); 
    }
    // Particle identification failed. Checking why...
    else {
      // Particle was identified as more than one particle type. 
      if(particleIdentified > 1) {
	G4cout <<"G4InclCascadeInterface: One outcoming particle was identified as";
	G4cout <<"more than one particle type. This is probably due to a bug in the interface." << G4endl;
	G4cout <<"Particle A:" << varntp->avv[particleI] << "Z: " << varntp->zvv[particleI] << G4endl;
	G4cout << "(particleIdentified =" << particleIdentified << ")"  << G4endl;
      }
    }
  }

  // End of conversion

  // Clean up: Clean up the number of generated particles in the
  // common block VARNTP_ for the processing of the next event.
  varntp->ntrack = 0;
  // End of cleanup.

// Free allocated memory
  delete varntp;
  delete abla;
  
  fillResult(cascadeParticles, theResult);
  return theResult;
}

void G4AblaEvaporation::fillResult( std::vector<G4DynamicParticle *> secondaryParticleVector,
				     G4FragmentVector * aResult )
{
  // Fill the vector pParticleChange with secondary particles stored in vector.
  G4cout <<"Size of the secondary particle vector = " << secondaryParticleVector.size() << G4endl;
  for ( size_t i = 0 ; i < secondaryParticleVector.size() ; i++ ) {
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
