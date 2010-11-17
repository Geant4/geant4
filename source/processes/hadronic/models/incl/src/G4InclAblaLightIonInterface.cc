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
// $Id: G4InclAblaLightIonInterface.cc,v 1.16 2010-11-17 20:19:09 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include <vector>

#include "G4InclAblaLightIonInterface.hh"
#include "G4FermiBreakUp.hh"
#include "math.h"
#include "G4GenericIon.hh"
#include "CLHEP/Random/Random.h"

G4InclAblaLightIonInterface::G4InclAblaLightIonInterface()
{
  hazard = new G4Hazard();

  const G4long* table_entry = CLHEP::HepRandom::getTheSeeds(); // Get random seed from CLHEP.
  hazard->ial = (*table_entry);

  varntp = new G4VarNtp();
  calincl = 0;
  ws = new G4Ws();
  mat = new G4Mat();
  incl = new G4Incl(hazard, calincl, ws, mat, varntp);
  useProjectileSpectator = true;
  useFermiBreakup = true;
  incl->setUseProjectileSpectators(useProjectileSpectator);
  if(!getenv("G4INCLABLANOFERMIBREAKUP")) { // Use Fermi Break-up by default if it is NOT explicitly disabled
    incl->setUseFermiBreakUp(true);
    useFermiBreakup = true;
  }
  verboseLevel = 0;
  if(getenv("G4INCLVERBOSE")) {
    verboseLevel = 1;
  }
}

G4InclAblaLightIonInterface::~G4InclAblaLightIonInterface()
{
  delete hazard;
  delete varntp;
  delete calincl;
  delete ws;
  delete mat;
  delete incl;
}

G4HadFinalState* G4InclAblaLightIonInterface::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  //  const G4bool useFermiBreakup = false;
  G4int maxTries = 200;

  G4int particleI;

  G4int baryonNumberBalanceInINCL = 0;
  G4int chargeNumberBalanceInINCL = 0;

  G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();

  // Increase the event number:
  eventNumber++;

  // Clean up the INCL input
  if(calincl != 0) {
    delete calincl;
    calincl = 0;
  }

  if (verboseLevel > 1) {
    G4cout << " >>> G4InclAblaLightIonInterface::ApplyYourself called" << G4endl;
  }

  if(verboseLevel > 1) {
    G4cout <<"G4InclAblaLightIonInterface: Now processing INCL4 event number:" << eventNumber << G4endl;
  }

  // Inverse kinematics for targets with Z = 1 and A = 1
  //  if(false) {
  G4LorentzRotation toBreit = aTrack.Get4Momentum().boostVector();

  if(theNucleus.GetZ_asInt() == 1 && theNucleus.GetA_asInt() == 1 && G4InclInput::canUseInverseKinematics(aTrack, theNucleus)) {
    G4ParticleDefinition *oldTargetDef = theTableOfParticles->GetIon(theNucleus.GetA_asInt(), theNucleus.GetZ_asInt(), 0.0);
    const G4ParticleDefinition *oldProjectileDef = aTrack.GetDefinition();

    if(oldProjectileDef != 0 && oldTargetDef != 0) {
    G4int oldTargetA = oldTargetDef->GetAtomicMass();
    G4int newTargetA = oldProjectileDef->GetAtomicMass();
    G4int newTargetZ = oldProjectileDef->GetAtomicNumber();

    if(newTargetA > 0 && newTargetZ > 0) {
      G4Nucleus swappedTarget(oldProjectileDef->GetAtomicMass(), oldProjectileDef->GetAtomicNumber());

      //      G4cout <<"Original projectile kinE = " << aTrack.GetKineticEnergy() / MeV << G4endl;

      // We need the same energy/nucleon.
      G4double projectileE = ((aTrack.GetKineticEnergy() / MeV) / newTargetA) * oldTargetA * MeV;

      //    G4cout <<"projectileE = " << projectileE << G4endl;
      G4DynamicParticle swappedProjectileParticle(oldTargetDef, G4ThreeVector(0.0, 0.0, 1.0), projectileE);
      const G4LorentzVector swapped4Momentum = (swappedProjectileParticle.Get4Momentum()*=toBreit);
      swappedProjectileParticle.Set4Momentum(swapped4Momentum);
      const G4HadProjectile swappedProjectile(swappedProjectileParticle);
      //  G4cout <<"New projectile kinE = " << swappedProjectile.GetKineticEnergy() / MeV << G4endl;
      calincl = new G4InclInput(swappedProjectile, swappedTarget, true);
    } else {
      G4cout <<"Badly defined target after swapping. Falling back to normal (non-swapped) mode." << G4endl;
      calincl = new G4InclInput(aTrack, theNucleus, false);
    }
    }
  } else {
    calincl = new G4InclInput(aTrack, theNucleus, false);
  }

  G4double eKin;
  G4double momx = 0.0, momy = 0.0, momz = 0.0;
  G4DynamicParticle *cascadeParticle = 0;
  G4ParticleDefinition *aParticleDefinition = 0;

  // INCL assumes the projectile particle is going in the direction of
  // the Z-axis. Here we construct proper rotation to convert the
  // momentum vectors of the outcoming particles to the original
  // coordinate system.
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();

  /*
  G4cout <<"Projectile theta = " << projectileMomentum.theta() << " phi = " << projectileMomentum.phi() << G4endl;
  G4cout <<"Projectile momentum "
	 << "(px = " << projectileMomentum.px()
	 << ", py = " << projectileMomentum.py()
	 << ", pz = " << projectileMomentum.pz() << ")" << G4endl;
  G4cout << "Projectile energy = " << bulletE << " MeV" << G4endl;
  */

  G4FermiBreakUp *fermiBreakUp = new G4FermiBreakUp();
  G4FragmentVector *theSpectatorFermiBreakupResult = 0;
  G4FragmentVector *theFermiBreakupResult = 0;

  theResult.Clear(); // Make sure the output data structure is clean.

  std::vector<G4DynamicParticle*> result; // Temporary list for the results

  // Map Geant4 particle types to corresponding INCL4 types.
  enum bulletParticleType {nucleus = 0, proton = 1, neutron = 2, pionPlus = 3, pionZero = 4, 
                           pionMinus = 5, deuteron = 6, triton = 7, he3 = 8, he4 = 9,
                           c12 = -12}; // Carbon beam support.

  G4int bulletType = calincl->bulletType();
  chargeNumberBalanceInINCL = calincl->targetZ();
  baryonNumberBalanceInINCL = calincl->targetA();

  //  G4cout <<"Type of the projectile (INCL projectile code): " << bulletType << G4endl;

  if(bulletType == proton) {
    chargeNumberBalanceInINCL += 1;
    baryonNumberBalanceInINCL += 1;
  } else if(bulletType == neutron) {
    baryonNumberBalanceInINCL += 1;
  } else if(bulletType == pionPlus) { //Note: positive pion doesn't contribute to the baryon and charge number counters
    chargeNumberBalanceInINCL += 1;
  } else if(bulletType == pionMinus) {
    chargeNumberBalanceInINCL -= 1;
  } else if(bulletType == deuteron) {
    chargeNumberBalanceInINCL += 1;
    baryonNumberBalanceInINCL += 2;
  } else if(bulletType == triton) {
    chargeNumberBalanceInINCL += 1;
    baryonNumberBalanceInINCL += 3;
  } else if(bulletType == he3) {
    chargeNumberBalanceInINCL += 2;
    baryonNumberBalanceInINCL += 3;
  } else if(bulletType == he4) {
    chargeNumberBalanceInINCL += 2;
    baryonNumberBalanceInINCL += 4;
  } if(bulletType == c12) {
    chargeNumberBalanceInINCL += 6;
    baryonNumberBalanceInINCL += 12;
  } if(bulletType == -666) {
    chargeNumberBalanceInINCL += calincl->extendedProjectileZ();
    baryonNumberBalanceInINCL += calincl->extendedProjectileA();
  }

  // Check wheter the input is acceptable.
  if((bulletType != 0) && ((calincl->targetA() != 1) && (calincl->targetZ() != 1))) {
    ws->nosurf = -2;  // Nucleus surface, -2 = Woods-Saxon 
    ws->xfoisa = 8;
    ws->npaulstr = 0;

    int nTries = 0;
    varntp->ntrack = 0;

    mat->nbmat = 1;
    mat->amat[0] = int(calincl->targetA());
    mat->zmat[0] = int(calincl->targetA());

    incl->setInput(calincl);
    incl->initIncl(true);

    while((varntp->ntrack <= 0) && (nTries < maxTries)) { // Loop until we produce real cascade
      nTries++;
      if(verboseLevel > 1) {
        G4cout <<"G4InclAblaLightIonInterface: Try number = " << nTries << G4endl; 
      }
      incl->processEventInclAbla(calincl, eventNumber);

      if(verboseLevel > 1) {
        G4cout <<"G4InclAblaLightIonInterface: number of tracks = " <<  varntp->ntrack <<G4endl;
      }
    }

    if(verboseLevel > 1) {
      /**
       * Diagnostic output
       */
      G4cout <<"G4InclAblaLightIonInterface: Bullet type: " << calincl->bulletType() << G4endl;
      G4cout <<"G4Incl4AblaCascadeInterface: Bullet energy: " << calincl->bulletE() << " MeV" << G4endl;
      if(bulletType == -666) {
	G4cout <<"   Extended projectile: A = " << calincl->extendedProjectileA()
	       <<" Z = " << calincl->extendedProjectileZ() << G4endl;
      }

      G4cout <<"G4InclAblaLightIonInterface: Target A:  " << calincl->targetA() << G4endl;
      G4cout <<"G4InclAblaLightIonInterface: Target Z:  " << calincl->targetZ() << G4endl;

      if(verboseLevel > 3) {
        diagdata <<"G4InclAblaLightIonInterface: Bullet type: " << calincl->bulletType() << G4endl;
        diagdata <<"G4InclAblaLightIonInterface: Bullet energy: " << calincl->bulletE() << " MeV" << G4endl;
        
        diagdata <<"G4InclAblaLightIonInterface: Target A:  " << calincl->targetA() << G4endl;
        diagdata <<"G4InclAblaLightIonInterface: Target Z:  " << calincl->targetZ() << G4endl;
      }
    }

    // Check whether a valid cascade was produced.
    // If not return the original bullet particle with the same momentum.
    if(varntp->ntrack <= 0) {
      if(verboseLevel > 1) {
        G4cout <<"WARNING G4InclAblaLightIonInterface: No cascade. Returning original particle with original momentum." << G4endl;
	G4cout <<"\t Reached maximum trials of 200 to produce inelastic scattering." << G4endl;
      }

      theResult.SetStatusChange(stopAndKill);
      
      if(bulletType == proton) {
        aParticleDefinition = G4Proton::ProtonDefinition();
      } else if(bulletType == neutron) {
        aParticleDefinition = G4Neutron::NeutronDefinition();
      } else if(bulletType == pionPlus) {
        aParticleDefinition = G4PionPlus::PionPlusDefinition();
      } else if(bulletType == pionZero) {
        aParticleDefinition = G4PionZero::PionZeroDefinition();
      } else if(bulletType == pionMinus) {
        aParticleDefinition = G4PionMinus::PionMinusDefinition();
      } else if(bulletType == deuteron) {
        aParticleDefinition = G4Deuteron::DeuteronDefinition();
      } else if(bulletType == triton) {
        aParticleDefinition = G4Triton::TritonDefinition();
      } else if(bulletType == he3) {
        aParticleDefinition = G4He3::He3Definition();
      } else if(bulletType == he4) {
        aParticleDefinition = G4Alpha::AlphaDefinition();
      } else { // Particle was not recognized. Probably an unsupported particle was given as input
	aParticleDefinition = 0;
      }

      if(aParticleDefinition != 0) {
	cascadeParticle = new G4DynamicParticle();
	cascadeParticle->SetDefinition(aParticleDefinition);
	cascadeParticle->Set4Momentum(aTrack.Get4Momentum());
	result.push_back(cascadeParticle);
      }
    }

    // Convert INCL4 output to Geant4 compatible data structures.
    // Elementary particles are converted to G4DynamicParticle.
    theResult.SetStatusChange(stopAndKill);
    
    for(particleI = 0; particleI <= varntp->ntrack; particleI++) { // Loop through the INCL4+ABLA output.
      // Get energy/momentum and construct momentum vector in INCL4 coordinates.
      //      if(varntp->itypcasc[particleI] == -1) continue; // Avoid nucleons that are part of the spectator
      if(varntp->avv[particleI] == 0 && varntp->zvv[particleI] == 0) continue;
      momx = varntp->plab[particleI]*std::sin(varntp->tetlab[particleI]*CLHEP::pi/180.0)*std::cos(varntp->philab[particleI]*CLHEP::pi/180.0)*MeV;
      momy = varntp->plab[particleI]*std::sin(varntp->tetlab[particleI]*CLHEP::pi/180.0)*std::sin(varntp->philab[particleI]*CLHEP::pi/180.0)*MeV;
      momz = varntp->plab[particleI]*std::cos(varntp->tetlab[particleI]*CLHEP::pi/180.0)*MeV;

      eKin = varntp->enerj[particleI] * MeV;

      G4ThreeVector momDirection(momx, momy, momz); // Direction of the particle.
      momDirection = momDirection.unit();
      if(verboseLevel > 2) {
	G4cout <<"G4InclAblaLightIonInterface: " << G4endl;
	G4cout <<"A    = " << varntp->avv[particleI] << " Z = "  << varntp->zvv[particleI] << G4endl;
	G4cout <<"eKin = " << eKin                   << " MeV"   << G4endl;
	G4cout <<"px   = " << momDirection.x()       << " py = " << momDirection.y()       <<" pz = " << momDirection.z() << G4endl;
      }

      G4int particleIdentified = 0; // Check particle ID.

      if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 1)) { // Proton
        cascadeParticle = 
          new G4DynamicParticle(G4Proton::ProtonDefinition(), momDirection, eKin);
        particleIdentified++;
	baryonNumberBalanceInINCL -= 1;
	chargeNumberBalanceInINCL -= 1;
      }

      if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 0)) { // Neutron
        cascadeParticle = 
          new G4DynamicParticle(G4Neutron::NeutronDefinition(), momDirection, eKin);
        particleIdentified++;
	baryonNumberBalanceInINCL -= 1;
      }

      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 1)) { // PionPlus
        cascadeParticle = 
          new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), momDirection, eKin);
        particleIdentified++;
	chargeNumberBalanceInINCL -= 1;
      }

      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 0)) { // PionZero
        cascadeParticle = 
          new G4DynamicParticle(G4PionZero::PionZeroDefinition(), momDirection, eKin);
        particleIdentified++;
	chargeNumberBalanceInINCL -= 0;
      }

      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == -1)) { // PionMinus
        cascadeParticle = 
          new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), momDirection, eKin);
        particleIdentified++;
	chargeNumberBalanceInINCL -= -1;
      }

      if((varntp->avv[particleI] > 1) && (varntp->zvv[particleI] >= 1)) { // Nucleus fragment
	G4ParticleDefinition * aIonDef = 0;

        G4int A = G4int(varntp->avv[particleI]);
        G4int Z = G4int(varntp->zvv[particleI]);
	G4double excitationE = G4double(varntp->exini) * MeV;

	if(verboseLevel > 1) {
	  G4cout <<"Finding ion: A = " << A << " Z = " << Z << " E* = " << excitationE/MeV << G4endl;
	}
	aIonDef = theTableOfParticles->GetIon(Z, A, excitationE);
	
	if(aIonDef == 0) {
	  if(verboseLevel > 1) {
	    G4cout <<"G4InclAblaLightIonInterface: " << G4endl;
	    G4cout <<"FATAL ERROR: aIonDef = 0" << G4endl;
	    G4cout <<"A = " << A << " Z = " << Z << " E* = " << excitationE << G4endl;
	  }
	}

	if(aIonDef != 0) { // If the ion was identified add it to output.
	  cascadeParticle =
	    new G4DynamicParticle(aIonDef, momDirection, eKin);
	  particleIdentified++;
	  baryonNumberBalanceInINCL -= A;
	  chargeNumberBalanceInINCL -= Z;
	}
      }
	
      if(particleIdentified == 1) { // Particle identified properly.
	cascadeParticle->Set4Momentum(cascadeParticle->Get4Momentum()*=toLabFrame);
        result.push_back(cascadeParticle);
      }
      else { // Particle identification failed.
        if(particleIdentified > 1) { // Particle was identified as more than one particle type. 
	  if(verboseLevel > 1) {
	    G4cout <<"G4InclAblaLightIonInterface: One outcoming particle was identified as";
	    G4cout <<"more than one particle type. This is probably due to a bug in the interface." << G4endl;
	    G4cout <<"Particle A:" << varntp->avv[particleI] << "Z: " << varntp->zvv[particleI] << G4endl;
	    G4cout << "(particleIdentified =" << particleIdentified << ")"  << G4endl;
	  }
        }
      }
    }

    // Spectator nucleus Fermi break-up
    if(useFermiBreakup && useProjectileSpectator && varntp->masp > 1) {
      baryonNumberBalanceInINCL -= G4int(varntp->masp);
      G4double nuclearMass = G4NucleiProperties::GetNuclearMass(G4int(varntp->masp), G4int(varntp->mzsp)) + varntp->exsp * MeV;
      // Use momentum scaling to compensate for different masses in G4 and INCL:
      G4double momentumScaling = G4InclUtils::calculate4MomentumScaling(G4int(varntp->masp),
									G4int(varntp->mzsp),
									varntp->exsp,
									varntp->spectatorT,
									varntp->spectatorP1,
									varntp->spectatorP2,
									varntp->spectatorP3);
      G4LorentzVector p4(momentumScaling * varntp->spectatorP1 * MeV, momentumScaling * varntp->spectatorP2 * MeV,
			 momentumScaling * varntp->spectatorP3 * MeV,
			 varntp->spectatorT * MeV + nuclearMass);
      // Four-momentum, baryon number and charge balance:
      G4LorentzVector fourMomentumBalance = p4;
      G4int baryonNumberBalance = G4int(varntp->masp);
      chargeNumberBalanceInINCL -= G4int(varntp->mzsp);
      G4int chargeBalance = G4int(varntp->mzsp);

      G4LorentzRotation toFragmentZ;
      // Assume that Fermi breakup uses Z as the direction of the projectile
      toFragmentZ.rotateZ(-p4.theta());
      toFragmentZ.rotateY(-p4.phi());
      G4LorentzRotation toFragmentLab = toFragmentZ.inverse();
      //      p4 *= toFragmentZ;
      
      G4LorentzVector p4rest = p4;
      //      p4rest.boost(-p4.boostVector());
      if(verboseLevel > 0) {
	G4cout <<"Spectator nucleus:" << G4endl;
	G4cout <<"p4: " << G4endl;
	G4cout <<" px: " << p4.px() <<" py: " << p4.py() <<" pz: " << p4.pz() << G4endl;
	G4cout <<" E = " << p4.e() << G4endl;
	G4cout <<"p4rest: " << G4endl;
	G4cout <<" px: " << p4rest.px() <<" py: " << p4rest.py() <<" pz: " << p4rest.pz() << G4endl;
	G4cout <<" E = " << p4rest.e() << G4endl;
      }
      G4Fragment theSpectatorNucleus(G4int(varntp->masp), G4int(varntp->mzsp), p4rest);
      theSpectatorFermiBreakupResult = fermiBreakUp->BreakItUp(theSpectatorNucleus);
      if(theSpectatorFermiBreakupResult != 0) {
      G4FragmentVector::iterator fragment;
      for(fragment = theSpectatorFermiBreakupResult->begin(); fragment != theSpectatorFermiBreakupResult->end(); fragment++) {
	G4ParticleDefinition *theFragmentDefinition = 0;
	if((*fragment)->GetA_asInt() == 1 && (*fragment)->GetZ_asInt() == 0) { // Neutron
	  theFragmentDefinition = G4Neutron::NeutronDefinition();
	} else if ((*fragment)->GetA_asInt() == 1 && (*fragment)->GetZ_asInt() == 1) {
	  theFragmentDefinition = G4Proton::ProtonDefinition();
	} else {
	  theFragmentDefinition = theTableOfParticles->GetIon((*fragment)->GetZ_asInt(), (*fragment)->GetA_asInt(), (*fragment)->GetExcitationEnergy());
	}
	if(theFragmentDefinition != 0) {
	  G4DynamicParticle *theFragment = new G4DynamicParticle(theFragmentDefinition, (*fragment)->GetMomentum());
	  G4LorentzVector labMomentum = theFragment->Get4Momentum();
	  //	  labMomentum.boost(p4.boostVector());
	  //	  labMomentum *= toFragmentLab;
	  //	  labMomentum *= toLabFrame;
	  theFragment->Set4Momentum(labMomentum);
	  fourMomentumBalance -= theFragment->Get4Momentum();
	  baryonNumberBalance -= theFragmentDefinition->GetAtomicMass();
	  chargeBalance -= theFragmentDefinition->GetAtomicNumber();
	  if(verboseLevel > 0) {
	    G4cout <<"Resulting fragment: " << G4endl;
	    G4cout <<" kinetic energy = " << theFragment->GetKineticEnergy() / MeV << " MeV" << G4endl;
	    G4cout <<" momentum = " << theFragment->GetMomentum().mag() / MeV << " MeV" << G4endl;
	  }
	  result.push_back(theFragment);
	} else {
	  G4cout <<"G4InclAblaCascadeInterface: Error. Fragment produced by Fermi break-up does not exist." 
		 << G4endl;
	  G4cout <<"Resulting fragment: " << G4endl;
	  G4cout <<" Z = " << (*fragment)->GetZ_asInt() << G4endl;
	  G4cout <<" A = " << (*fragment)->GetA_asInt() << G4endl;
	  G4cout <<" Excitation : " << (*fragment)->GetExcitationEnergy() / MeV << " MeV" << G4endl;
	  G4cout <<" momentum = " << (*fragment)->GetMomentum().mag() / MeV << " MeV" << G4endl;
	}
      }
      delete theSpectatorFermiBreakupResult;
      theSpectatorFermiBreakupResult = 0;

      if(std::abs(fourMomentumBalance.mag() / MeV) > 0.1 * MeV) { 
	G4cout <<"Four-momentum balance after spectator nucleus Fermi break-up:" << G4endl;
	G4cout <<"Magnitude: " << fourMomentumBalance.mag() / MeV << " MeV" << G4endl;
	G4cout <<"Vector components (px, py, pz, E) = ("
	       << fourMomentumBalance.px() << ", "
	       << fourMomentumBalance.py() << ", "
	       << fourMomentumBalance.pz() << ", "
	       << fourMomentumBalance.e() << ")" << G4endl;
      }
      if(baryonNumberBalance != 0) {
	G4cout <<"Event " << eventNumber << ": Baryon number balance after spectator nucleus Fermi break-up: " << baryonNumberBalance << G4endl;
      }
      if(chargeBalance != 0) {
	G4cout <<"Event " << eventNumber <<": Charge balance after spectator nucleus Fermi break-up: " << chargeBalance << G4endl;
      }
    }
  }

    // Finally do Fermi break-up if needed
    if(varntp->needsFermiBreakup && varntp->massini > 0) {
      baryonNumberBalanceInINCL -= G4int(varntp->massini);
      chargeNumberBalanceInINCL -= G4int(varntp->mzini);
      // Call Fermi Break-up
      G4double nuclearMass = G4NucleiProperties::GetNuclearMass(G4int(varntp->massini), G4int(varntp->mzini)) + varntp->exini * MeV;
      G4LorentzVector fragmentMomentum(varntp->pxrem * MeV, varntp->pyrem * MeV, varntp->pzrem * MeV,
				       varntp->erecrem * MeV + nuclearMass);
      G4double momentumScaling = G4InclUtils::calculate4MomentumScaling(G4int(varntp->massini), G4int(varntp->mzini),
									varntp->exini,
									varntp->erecrem,
									varntp->pxrem,
									varntp->pyrem,
									varntp->pzrem);
      G4LorentzVector p4(momentumScaling * varntp->pxrem * MeV, momentumScaling * varntp->pyrem * MeV,
			 momentumScaling * varntp->pzrem * MeV,
      			 varntp->erecrem + nuclearMass);

      // For four-momentum, baryon number and charge conservation check:
      G4LorentzVector fourMomentumBalance = p4;
      G4int baryonNumberBalance = G4int(varntp->massini);
      G4int chargeBalance = G4int(varntp->mzini);

      G4LorentzRotation toFragmentZ;
      toFragmentZ.rotateZ(-p4.theta());
      toFragmentZ.rotateY(-p4.phi());
      G4LorentzRotation toFragmentLab = toFragmentZ.inverse();
      //      p4 *= toFragmentZ;

      G4LorentzVector p4rest = p4;
      //      p4rest.boost(-p4.boostVector());
      if(verboseLevel > 0) {
	G4cout <<"Cascade remnant nucleus:" << G4endl;
	G4cout <<"p4: " << G4endl;
	G4cout <<" px: " << p4.px() <<" py: " << p4.py() <<" pz: " << p4.pz() << G4endl;
	G4cout <<" E = " << p4.e() << G4endl;

	G4cout <<"p4rest: " << G4endl;
	G4cout <<" px: " << p4rest.px() <<" py: " << p4rest.py() <<" pz: " << p4rest.pz() << G4endl;
	G4cout <<" E = " << p4rest.e() << G4endl;
      }

      G4Fragment theCascadeRemnant(G4int(varntp->massini), G4int(varntp->mzini), p4rest);
      theFermiBreakupResult = fermiBreakUp->BreakItUp(theCascadeRemnant);
      if(theFermiBreakupResult != 0) {
      G4FragmentVector::iterator fragment;
      for(fragment = theFermiBreakupResult->begin(); fragment != theFermiBreakupResult->end(); fragment++) {
	G4ParticleDefinition *theFragmentDefinition = 0;
	if((*fragment)->GetA_asInt() == 1 && (*fragment)->GetZ_asInt() == 0) { // Neutron
	  theFragmentDefinition = G4Neutron::NeutronDefinition();
	} else if ((*fragment)->GetA_asInt() == 1 && (*fragment)->GetZ_asInt() == 1) {
	  theFragmentDefinition = G4Proton::ProtonDefinition();
	} else {
	  theFragmentDefinition = theTableOfParticles->GetIon((*fragment)->GetZ_asInt(), (*fragment)->GetA_asInt(), (*fragment)->GetExcitationEnergy());
	}

	if(theFragmentDefinition != 0) {
	  G4DynamicParticle *theFragment = new G4DynamicParticle(theFragmentDefinition, (*fragment)->GetMomentum());
	  G4LorentzVector labMomentum = theFragment->Get4Momentum();
	  //	  labMomentum.boost(p4.boostVector());
	  //	  labMomentum *= toFragmentLab;
	  //	  labMomentum *= toLabFrame;
	  theFragment->Set4Momentum(labMomentum);
	  fourMomentumBalance -= theFragment->Get4Momentum();
	  baryonNumberBalance -= theFragmentDefinition->GetAtomicMass();
	  chargeBalance -= theFragmentDefinition->GetAtomicNumber();
	  if(verboseLevel > 0) {
	    G4cout <<"Resulting fragment: " << G4endl;
	    G4cout <<" kinetic energy = " << theFragment->GetKineticEnergy() / MeV << " MeV" << G4endl;
	    G4cout <<" momentum = " << theFragment->GetMomentum().mag() / MeV << " MeV" << G4endl;
	  }
	  result.push_back(theFragment);
	} else {
	  G4cout <<"G4InclAblaCascadeInterface: Error. Fragment produced by Fermi break-up does not exist." << G4endl;
	  G4cout <<"Resulting fragment: " << G4endl;
	  G4cout <<" Z = " << (*fragment)->GetZ_asInt() << G4endl;
	  G4cout <<" A = " << (*fragment)->GetA_asInt() << G4endl;
	  G4cout <<" Excitation : " << (*fragment)->GetExcitationEnergy() / MeV << " MeV" << G4endl;
	  G4cout <<" momentum = " << (*fragment)->GetMomentum().mag() / MeV << " MeV" << G4endl;
	}
      }
      delete theFermiBreakupResult;
      theFermiBreakupResult = 0;

      if(std::abs(fourMomentumBalance.mag() / MeV) > 0.1 * MeV) { 
	G4cout <<"Four-momentum balance after remnant nucleus Fermi break-up:" << G4endl;
	G4cout <<"Magnitude: " << fourMomentumBalance.mag() / MeV << " MeV" << G4endl;
	G4cout <<"Vector components (px, py, pz, E) = ("
	       << fourMomentumBalance.px() << ", "
	       << fourMomentumBalance.py() << ", "
	       << fourMomentumBalance.pz() << ", "
	       << fourMomentumBalance.e() << ")" << G4endl;
      }
      if(baryonNumberBalance != 0) {
	G4cout <<"Baryon number balance after remnant nucleus Fermi break-up: " << baryonNumberBalance << G4endl;
      }
      if(chargeBalance != 0) {
	G4cout <<"Charge balance after remnant nucleus Fermi break-up: " << chargeBalance << G4endl;
      }
    }
    }

    varntp->ntrack = 0; // Clean up the number of generated particles in the event.

    if(baryonNumberBalanceInINCL != 0 && verboseLevel > 1) {
      G4cout <<"Event " << eventNumber <<": G4InclAblaLightIonInterface: Baryon number conservation problem in INCL detected!" << G4endl;
      G4cout <<"Baryon number balance: " << baryonNumberBalanceInINCL << G4endl;
      if(baryonNumberBalanceInINCL < 0) {
	G4cout <<"Event " << eventNumber <<": Too many outcoming baryons!" << G4endl;
      } else if(baryonNumberBalanceInINCL > 0) {
	G4cout <<"Event " << eventNumber <<": Too few outcoming baryons!" << G4endl;
      }
    }

    if(chargeNumberBalanceInINCL != 0 && verboseLevel > 1) {
      G4cout <<"Event " << eventNumber <<": G4InclAblaLightIonInterface: Charge number conservation problem in INCL detected!" << G4endl;
      G4cout <<"Event " << eventNumber <<": Charge number balance: " << chargeNumberBalanceInINCL << G4endl;
    }
  }
  /**
   * Report unsupported features.
   * (Check bullet, target, energy range)
   */
  else { // If the bullet type was not recognized by the interface, it will be returned back without any interaction.
    theResult.SetStatusChange(stopAndKill);

    G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();
    cascadeParticle = new G4DynamicParticle(theTableOfParticles->FindParticle(aTrack.GetDefinition()), aTrack.Get4Momentum());

    result.push_back(cascadeParticle);

    if(verboseLevel > 1) {
      G4cout <<"G4InclAblaLightIonInterface: Error processing event number (internal) " << eventNumber << G4endl;
    }
    if(verboseLevel > 3) {
      diagdata <<"G4InclAblaLightIonInterface: Error processing event number (internal) " << eventNumber << G4endl;
    }

    if(bulletType == 0) {
      if(verboseLevel > 1) {
	G4cout <<"G4InclAblaLightIonInterface: Unknown bullet type" << G4endl;      
	G4cout <<"Bullet particle name: " << cascadeParticle->GetDefinition()->GetParticleName() << G4endl;
      }
      if(verboseLevel > 3) {
        diagdata <<"G4InclAblaLightIonInterface: Unknown bullet type" << G4endl;      
        diagdata <<"Bullet particle name: " << cascadeParticle->GetDefinition()->GetParticleName() << G4endl;
      }
    }

    if((calincl->targetA() == 1) && (calincl->targetZ() == 1)) { // Unsupported target
      if(verboseLevel > 1) {
	G4cout <<"Unsupported target: " << G4endl;
	G4cout <<"Target A: " << calincl->targetA() << G4endl;
	G4cout <<"TargetZ: " << calincl->targetZ() << G4endl;
      }
      if(verboseLevel > 3) {
        diagdata <<"Unsupported target: " << G4endl;
        diagdata <<"Target A: " << calincl->targetA() << G4endl;
	diagdata <<"TargetZ: " << calincl->targetZ() << G4endl;
      }
    }

    if(calincl->bulletE() < 100) { // INCL does not support E < 100 MeV.
      if(verboseLevel > 1) {
	G4cout <<"Unsupported bullet energy: " << calincl->bulletE() << " MeV. (Lower limit is 100 MeV)." << G4endl;
	G4cout <<"WARNING: Returning the original bullet with original energy back to Geant4." << G4endl;
      }
      if(verboseLevel > 3) {
        diagdata <<"Unsupported bullet energy: " << calincl->bulletE() << " MeV. (Lower limit is 100 MeV)." << G4endl;
      }
    }

    if(verboseLevel > 3) {
      diagdata <<"WARNING: returning the original bullet with original energy back to Geant4." << G4endl;
    }
  }

  // Finally copy the accumulated secondaries into the result collection:
  G4ThreeVector boostVector = aTrack.Get4Momentum().boostVector();
  G4LorentzRotation boostBack = toBreit.inverse();

  for(std::vector<G4DynamicParticle*>::iterator i = result.begin(); i != result.end(); ++i) {
    // If the calculation was performed in inverse kinematics we have to
    // convert the result back...
    if(calincl->isInverseKinematics()) {
      G4LorentzVector mom = (*i)->Get4Momentum();
      mom.setPz(-1.0 * mom.pz()); // Reverse the z-component of the momentum vector
      mom *= boostBack;
      (*i)->Set4Momentum(mom);
    }
    theResult.AddSecondary((*i));
  }

  delete fermiBreakUp;
  delete calincl;
  calincl = 0;
  return &theResult;
} 

G4ReactionProductVector* G4InclAblaLightIonInterface::Propagate(G4KineticTrackVector* , G4V3DNucleus* ) {
  return 0;
}


