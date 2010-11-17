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
// $Id: G4InclCascadeInterface.cc,v 1.15 2010-11-17 20:19:09 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

//#define DEBUGINCL 1

#include "G4InclCascadeInterface.hh"
#include "G4FermiBreakUp.hh"
#include "math.h"
#include "G4GenericIon.hh"
#include "CLHEP/Random/Random.h"

G4InclCascadeInterface::G4InclCascadeInterface(const G4String& nam)
  :G4VIntraNuclearTransportModel(nam)
{
  hazard = new G4Hazard();
  const G4long* table_entry = CLHEP::HepRandom::getTheSeeds(); // Get random seed from CLHEP.
  hazard->ial = (*table_entry);

  varntp = new G4VarNtp();
  calincl = 0;
  ws = new G4Ws();
  mat = new G4Mat();
  incl = new G4Incl(hazard, calincl, ws, mat, varntp);

  theExcitationHandler = new G4ExcitationHandler;
  thePrecoModel = new G4PreCompoundModel(theExcitationHandler);

  if(!getenv("G4INCLABLANOFERMIBREAKUP")) { // Use Fermi Break-up by default if it is NOT explicitly disabled
    incl->setUseFermiBreakUp(true);
  }

  verboseLevel = 0;
}

G4InclCascadeInterface::~G4InclCascadeInterface()
{
  delete thePrecoModel;
  delete theExcitationHandler;

  delete hazard;
  delete varntp;
  delete ws;
  delete mat;
  delete incl;
}

G4HadFinalState* G4InclCascadeInterface::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  G4int maxTries = 200;

  G4int particleI;
  G4int bulletType = 0;

  // Print diagnostic messages: 0 = silent, 1 and 2 = verbose
  verboseLevel = 0;

  // Increase the event number:
  eventNumber++;

  if (verboseLevel > 1) {
    G4cout << " >>> G4InclCascadeInterface::ApplyYourself called" << G4endl;
  }

  if(verboseLevel > 1) {
    G4cout <<"G4InclCascadeInterface: Now processing INCL4 event number:" << eventNumber << G4endl;
  }

#ifdef DEBUGINCL
  G4cout <<"Bullet energy = " << bulletE / MeV << G4endl;
#endif

  G4double eKin;
  G4double momx = 0.0, momy = 0.0, momz = 0.0;
  G4DynamicParticle *cascadeParticle = 0;
  G4ParticleDefinition *aParticleDefinition = 0;
  
  G4ReactionProductVector *thePrecoResult = 0;
  G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();

  // INCL assumes the projectile particle is going in the direction of
  // the Z-axis. Here we construct proper rotation to convert the
  // momentum vectors of the outcoming particles to the original
  // coordinate system.
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();

  theResult.Clear(); // Make sure the output data structure is clean.

  calincl = new G4InclInput(aTrack, theNucleus, false);
  incl->setInput(calincl);

  //  G4InclInput::printProjectileTargetInfo(aTrack, theNucleus);
  //  calincl->printInfo();

#ifdef DEBUGINCL
  G4int baryonBullet = 0, chargeBullet = 0;
  if(bulletType == proton || bulletType == neutron) baryonBullet = 1;
  if(bulletType == proton || bulletType == pionPlus) chargeBullet = 1;
  if(bulletType == pionMinus) chargeBullet = -1;
  G4int baryonNumber = int(std::floor(targetA)) + baryonBullet;
  G4int chargeNumber = int(std::floor(targetZ)) + chargeBullet;  
  G4double mass = aTrack.GetDefinition()->GetPDGMass();
  G4double amass = theNucleus.AtomicMass(targetA, targetZ);
  G4double eKinSum = bulletE;
  G4LorentzVector labv = G4LorentzVector(0.0, 0.0, std::sqrt(bulletE*(bulletE + 2.*mass)), bulletE + mass + amass);
  G4LorentzVector labvA = G4LorentzVector(0.0, 0.0, 0.0, 0.0);
  G4cout <<"Energy in the beginning = " << labv.e() / MeV << G4endl;
#endif

  // Check wheter the input is acceptable.
  if((calincl->bulletType() != 0) && ((calincl->targetA() != 1) && (calincl->targetZ() != 1))) {
    ws->nosurf = -2;  // Nucleus surface, -2 = Woods-Saxon 
    ws->xfoisa = 8;
    ws->npaulstr = 0;

    int nTries = 0;
    varntp->ntrack = 0;

    mat->nbmat = 1;
    mat->amat[0] = int(calincl->targetA());
    mat->zmat[0] = int(calincl->targetZ());

    incl->initIncl(true);

    while((varntp->ntrack <= 0) && (nTries < maxTries)) { // Loop until we produce real cascade
      nTries++;
      if(verboseLevel > 1) {
        G4cout <<"G4InclCascadeInterface: Try number = " << nTries << G4endl; 
      }
      incl->processEventIncl(calincl);

      if(verboseLevel > 1) {
        G4cout <<"G4InclCascadeInterface: number of tracks = " <<  varntp->ntrack <<G4endl;
      }
    }

    if(verboseLevel > 1) {
      /**
       * Diagnostic output
       */
      G4cout <<"G4InclCascadeInterface: Bullet type: " << calincl->bulletType() << G4endl;
      G4cout <<"G4Incl4AblaCascadeInterface: Bullet energy: " << calincl->bulletE() << " MeV" << G4endl;

      G4cout <<"G4InclCascadeInterface: Target A:  " << calincl->targetA() << G4endl;
      G4cout <<"G4InclCascadeInterface: Target Z:  " << calincl->targetZ() << G4endl;

      if(verboseLevel > 3) {
        diagdata <<"G4InclCascadeInterface: Bullet type: " << calincl->bulletType() << G4endl;
        diagdata <<"G4InclCascadeInterface: Bullet energy: " << calincl->bulletE() << " MeV" << G4endl;
        
        diagdata <<"G4InclCascadeInterface: Target A:  " << calincl->targetA() << G4endl;
        diagdata <<"G4InclCascadeInterface: Target Z:  " << calincl->targetZ() << G4endl;
      }

    }

    // Check whether a valid cascade was produced.
    // If not return the original bullet particle with the same momentum.
    if(varntp->ntrack <= 0) {
      if(verboseLevel > 1) {
        G4cout <<"WARNING G4InclCascadeInterface: No cascade. Returning original particle with original momentum." << G4endl;
	G4cout <<"\t Reached maximum trials of 200 to produce inelastic scattering." << G4endl;
      }

      theResult.SetStatusChange(stopAndKill);
      
      G4int bulletType = calincl->bulletType();
      aParticleDefinition = G4InclInput::getParticleDefinition(bulletType);

      if(aParticleDefinition != 0) {
	cascadeParticle = new G4DynamicParticle();
	cascadeParticle->SetDefinition(aParticleDefinition);
	cascadeParticle->Set4Momentum(aTrack.Get4Momentum());
	theResult.AddSecondary(cascadeParticle);
      }
    }

    // Convert INCL4 output to Geant4 compatible data structures.
    // Elementary particles are converted to G4DynamicParticle.
    theResult.SetStatusChange(stopAndKill);

#ifdef DEBUGINCL
    G4cout << "E [MeV]" << std::setw(12)
	   << " Ekin [MeV]" << std::setw(12)
	   << "Px [MeV]" << std::setw(12)
	   << " Py [MeV]" << std::setw(12)
	   << "Pz [MeV]" << std::setw(12)
	   << "Pt [MeV]" << std::setw(12)
	   << "A" << std::setw(12)
	   << "Z" << G4endl;
#endif

    for(particleI = 0; particleI < varntp->ntrack; particleI++) { // Loop through the INCL4+ABLA output.
      // Get energy/momentum and construct momentum vector in INCL4 coordinates.
      momx = varntp->plab[particleI]*std::sin(varntp->tetlab[particleI]*CLHEP::pi/180.0)*std::cos(varntp->philab[particleI]*CLHEP::pi/180.0)*MeV;
      momy = varntp->plab[particleI]*std::sin(varntp->tetlab[particleI]*CLHEP::pi/180.0)*std::sin(varntp->philab[particleI]*CLHEP::pi/180.0)*MeV;
      momz = varntp->plab[particleI]*std::cos(varntp->tetlab[particleI]*CLHEP::pi/180.0)*MeV;

      eKin = varntp->enerj[particleI] * MeV;

      G4ThreeVector momDirection(momx, momy, momz); // Direction of the particle.
      momDirection = momDirection.unit();
      if(verboseLevel > 2) {
	G4cout <<"G4InclCascadeInterface: " << G4endl;
	G4cout <<"A    = " << varntp->avv[particleI] << " Z = "  << varntp->zvv[particleI] << G4endl;
	G4cout <<"eKin = " << eKin                   << " MeV"   << G4endl;
	G4cout <<"px   = " << momDirection.x()       << " py = " << momDirection.y()       <<" pz = " << momDirection.z() << G4endl;
      }

      G4int particleIdentified = 0; // Check particle ID.

      if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 1)) { // Proton
        cascadeParticle = 
          new G4DynamicParticle(G4Proton::ProtonDefinition(), momDirection, eKin);
        particleIdentified++;
#ifdef DEBUGINCL
	baryonNumber--;
	chargeNumber--;
#endif
      }

      if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 0)) { // Neutron
        cascadeParticle = 
          new G4DynamicParticle(G4Neutron::NeutronDefinition(), momDirection, eKin);
        particleIdentified++;
#ifdef DEBUGINCL
	baryonNumber--;
#endif
      }

      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 1)) { // PionPlus
        cascadeParticle = 
          new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), momDirection, eKin);
        particleIdentified++;
#ifdef DEBUGINCL
	chargeNumber--;
#endif
      }

      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 0)) { // PionZero
        cascadeParticle = 
          new G4DynamicParticle(G4PionZero::PionZeroDefinition(), momDirection, eKin);
        particleIdentified++;
      }

      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == -1)) { // PionMinus
        cascadeParticle = 
          new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), momDirection, eKin);
        particleIdentified++;
#ifdef DEBUGINCL
	chargeNumber++;
#endif
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
	    G4cout <<"G4InclCascadeInterface: " << G4endl;
	    G4cout <<"FATAL ERROR: aIonDef = 0" << G4endl;
	    G4cout <<"A = " << A << " Z = " << Z << " E* = " << excitationE << G4endl;
	  }
	}

	if(aIonDef != 0) { // If the ion was identified add it to output.
	  cascadeParticle =
	    new G4DynamicParticle(aIonDef, momDirection, eKin);
	  particleIdentified++;
#ifdef DEBUGINCL
	  baryonNumber = baryonNumber - A;
	  chargeNumber = chargeNumber - Z;
#endif
	}
      }
	
      if(particleIdentified == 1) { // Particle identified properly.
        cascadeParticle->Set4Momentum(cascadeParticle->Get4Momentum()*=toLabFrame);
#ifdef DEBUGINCL
	G4ParticleDefinition *pd  = cascadeParticle->GetDefinition();
	G4LorentzVector fm  = cascadeParticle->Get4Momentum();
	G4ThreeVector mom = cascadeParticle->GetMomentum();
	G4double m = pd->GetPDGMass();
	G4double p = mom.mag();
	labv -= fm;
	if(varntp->avv[particleI] > 1) {
	  labvA += fm;
	}
	G4double px = mom.x() * MeV;
	G4double py = mom.y() * MeV;
	G4double pz = mom.z() * MeV;
	G4double pt = std::sqrt(px*px+py*py);
	G4double e  = fm.e();
	eKinSum -= cascadeParticle->GetKineticEnergy() * MeV;
	G4double exE;
	if(varntp->avv[particleI] > 1) {
	  exE = varntp->exini;
	}
	else {
	  exE = 0.0;
	}
	G4cout << fm.e() / MeV
	       << std::setw(12) << cascadeParticle->GetKineticEnergy() / MeV
	       << std::setw(12) << mom.x() / MeV
	       << std::setw(12) << mom.y() / MeV
	       << std::setw(12) << mom.z() / MeV
	       << std::setw(12) << pt / MeV
	       << std::setw(12) << varntp->avv[particleI]
	       << std::setw(12) << varntp->zvv[particleI] << G4endl;
#endif
        theResult.AddSecondary(cascadeParticle); // Put data into G4HadFinalState.
      }
      else { // Particle identification failed.
        if(particleIdentified > 1) { // Particle was identified as more than one particle type. 
	  if(verboseLevel > 1) {
	    G4cout <<"G4InclCascadeInterface: One outcoming particle was identified as";
	    G4cout <<"more than one particle type. This is probably due to a bug in the interface." << G4endl;
	    G4cout <<"Particle A:" << varntp->avv[particleI] << "Z: " << varntp->zvv[particleI] << G4endl;
	    G4cout << "(particleIdentified =" << particleIdentified << ")"  << G4endl;
	  }
        }
      }
    }

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
      p4 *= toFragmentZ;

      G4LorentzVector p4rest = p4;
      p4rest.boost(-p4.boostVector());
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
      thePrecoResult = thePrecoModel->DeExcite(theCascadeRemnant);
      if(thePrecoResult != 0) {
      G4ReactionProductVector::iterator fragment;
      for(fragment = thePrecoResult->begin(); fragment != thePrecoResult->end(); fragment++) {
	G4ParticleDefinition *theFragmentDefinition = (*fragment)->GetDefinition();

	if(theFragmentDefinition != 0) {
	  G4DynamicParticle *theFragment = new G4DynamicParticle(theFragmentDefinition, (*fragment)->GetMomentum());
	  G4LorentzVector labMomentum = theFragment->Get4Momentum();
	  labMomentum.boost(p4.boostVector());
	  labMomentum *= toFragmentLab;
	  labMomentum *= toLabFrame;
	  theFragment->Set4Momentum(labMomentum);
	  fourMomentumBalance -= theFragment->Get4Momentum();
	  baryonNumberBalance -= theFragmentDefinition->GetAtomicMass();
	  chargeBalance -= theFragmentDefinition->GetAtomicNumber();
	  if(verboseLevel > 0) {
	    G4cout <<"Resulting fragment: " << G4endl;
	    G4cout <<" kinetic energy = " << theFragment->GetKineticEnergy() / MeV << " MeV" << G4endl;
	    G4cout <<" momentum = " << theFragment->GetMomentum().mag() / MeV << " MeV" << G4endl;
	  }
	  theResult.AddSecondary(theFragment);
	} else {
	  G4cout <<"G4InclCascadeInterface: Error. Fragment produced by Fermi break-up does not exist." << G4endl;
	  G4cout <<"Resulting fragment: " << G4endl;
	  G4cout <<" momentum = " << (*fragment)->GetMomentum().mag() / MeV << " MeV" << G4endl;
	}
      }
      delete thePrecoResult;
      thePrecoResult = 0;

      if(verboseLevel > 1 && std::abs(fourMomentumBalance.mag() / MeV) > 0.1 * MeV) { 
	G4cout <<"Four-momentum balance after remnant nucleus Fermi break-up:" << G4endl;
	G4cout <<"Magnitude: " << fourMomentumBalance.mag() / MeV << " MeV" << G4endl;
	G4cout <<"Vector components (px, py, pz, E) = ("
	       << fourMomentumBalance.px() << ", "
	       << fourMomentumBalance.py() << ", "
	       << fourMomentumBalance.pz() << ", "
	       << fourMomentumBalance.e() << ")" << G4endl;
      }
      if(baryonNumberBalance != 0 && verboseLevel > 1) {
	G4cout <<"Baryon number balance after remnant nucleus Fermi break-up: " << baryonNumberBalance << G4endl;
      }
      if(chargeBalance != 0 && verboseLevel > 1) {
	G4cout <<"Charge balance after remnant nucleus Fermi break-up: " << chargeBalance << G4endl;
      }
      }
      //    } // if(needsFermiBreakUp)

#ifdef DEBUGINCL
    G4cout <<"--------------------------------------------------------------------------------" << G4endl;
    G4double pt = std::sqrt(std::pow(labv.x(), 2) + std::pow(labv.y(), 2));
    G4double ptA = std::sqrt(std::pow(labvA.x(), 2) + std::pow(labvA.y(), 2));
    G4cout << labv.e() / MeV << std::setw(12)
	   << eKinSum  / MeV << std::setw(12)
	   << labv.x() / MeV << std::setw(12)
	   << labv.y() / MeV << std::setw(12)
	   << labv.z() / MeV << std::setw(12)
	   << pt       / MeV << std::setw(12)
	   << baryonNumber << std::setw(12)
	   << chargeNumber << " totals" << G4endl;
    G4cout << " - " << std::setw(12)
	   << " - " << std::setw(12)
	   << labvA.x() / MeV << std::setw(12)
	   << labvA.y() / MeV << std::setw(12)
	   << labvA.z() / MeV << std::setw(12)
	   << ptA       / MeV << std::setw(12)
	   << " - " << std::setw(12) << " - " << " totals ABLA" << G4endl;
    G4cout << G4endl;
 
    if(verboseLevel > 3) {
      if(baryonNumber != 0) {
	G4cout <<"WARNING G4InclCascadeInterface: Baryon number conservation violated." << G4endl;
	G4cout <<"Baryon number balance after the event: " << baryonNumber << G4endl;
	if(baryonNumber < 0) {
	  G4cout <<"Too many baryons produced." << G4endl;
	} else {
	  G4cout <<"Too few baryons produced." << G4endl;
	}
      }
    }
#endif
   
    varntp->ntrack = 0; // Clean up the number of generated particles in the event.
  }
  /**
   * Report unsupported features.
   * (Check bullet, target, energy range)
   */
  else { // If the bullet type was not recognized by the interface, it will be returned back without any interaction.
    theResult.SetStatusChange(stopAndKill);

    cascadeParticle = new G4DynamicParticle(theTableOfParticles->FindParticle(aTrack.GetDefinition()), aTrack.Get4Momentum());

    theResult.AddSecondary(cascadeParticle);

    if(verboseLevel > 1) {
      G4cout <<"ERROR G4InclCascadeInterface: Processing event number (internal) failed " << eventNumber << G4endl;
    }
    if(verboseLevel > 3) {
      diagdata <<"ERROR G4InclCascadeInterface: Error processing event number (internal) failed " << eventNumber << G4endl;
    }

    if(bulletType == 0) {
      if(verboseLevel > 1) {
	G4cout <<"G4InclCascadeInterface: Unknown bullet type" << G4endl;      
	G4cout <<"Bullet particle name: " << cascadeParticle->GetDefinition()->GetParticleName() << G4endl;
      }
      if(verboseLevel > 3) {
        diagdata <<"G4InclCascadeInterface: Unknown bullet type" << G4endl;      
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

  delete calincl;
  calincl = 0;
  return &theResult;
} 

G4ReactionProductVector* G4InclCascadeInterface::Propagate(G4KineticTrackVector* , G4V3DNucleus* ) {
  return 0;
}


