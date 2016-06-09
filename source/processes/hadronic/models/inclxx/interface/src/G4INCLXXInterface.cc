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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLXXInterface.hh"
#include "G4INCLXXFactory.hh"
#include "math.h"
#include "G4GenericIon.hh"
#include "CLHEP/Random/Random.h"
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

G4INCLXXInterface::G4INCLXXInterface(const G4String& nam)
  :G4VIntraNuclearTransportModel(nam)
{
  // Use the environment variable G4INCLXX_NO_DE_EXCITATION to disable de-excitation
  if(getenv("G4INCLXX_NO_DE_EXCITATION")) {
    G4cout <<"INCL++ Interface: WARNING: De-excitation is completely disabled!" << G4endl;
    theExcitationHandler = 0;
  } else {
    theExcitationHandler = new G4ExcitationHandler;
  }

  if(getenv("G4INCLXX_STORE_RAW_DEBUG_OUTPUT")) {
    storeDebugOutput = true;
    debugOutputFile = new std::ofstream("inclDebug.out");
  } else {
    storeDebugOutput = false;
  }

  dumpInput = false;
}

G4INCLXXInterface::~G4INCLXXInterface()
{
  delete theExcitationHandler;
  if(storeDebugOutput) {
    debugOutputFile->close();
  }
}

G4HadFinalState* G4INCLXXInterface::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  const G4int maxTries = 200;

  //  G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();

  // Create new INCL config object that contains the projectile and target information. This object also contains the model parameters (e.g. 
  

  // INCL assumes the projectile particle is going in the direction of
  // the Z-axis. Here we construct proper rotation to convert the
  // momentum vectors of the outcoming particles to the original
  // coordinate system.
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum();

  // INCL++ assumes that the projectile is going in the direction of
  // the z-axis. In principle, if the coordinate system used by G4
  // hadronic framework is defined differently we need a rotation to
  // transform the INCL++ reaction products to the appropriate
  // frame. Please note that it isn't necessary to apply this
  // transform to the projectile because when creating the INCL++
  // projectile particle G4INCLXXFactory::createProjectile needs to
  // use only the projectile energy (direction is simply assumed to be
  // along z-axis).
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();
  // However, it turns out that the projectile given to us by G4
  // hadronic framework is already going in the direction of the
  // z-axis so this rotation is actually unnecessary. Both toZ and
  // toLabFrame turn out to be unit matrices as can be seen by
  // uncommenting the folowing two lines:
  // G4cout <<"toZ = " << toZ << G4endl;
  // G4cout <<"toLabFrame = " << toLabFrame << G4endl;

  theResult.Clear(); // Make sure the output data structure is clean.
  theResult.SetStatusChange(stopAndKill);

  std::list<G4Fragment> remnants;

  G4int nTries = 0;
  // INCL can generate transparent events. However, this is meaningful
  // only in the standalone code. In Geant4 we must "force" INCL to
  // produce a valid cascade.
  G4bool eventIsOK = false;
  do {
    G4INCL::INCL *theINCLModel = G4INCLXXFactory::createModel(theNucleus);
    G4INCL::Particle *theProjectile = G4INCLXXFactory::createProjectile(aTrack);

    if(dumpInput) {
      G4cout << theINCLModel->configToString() << G4endl;
    }
    const G4INCL::EventInfo eventInfo = theINCLModel->processEvent(theProjectile);
    //    eventIsOK = !eventInfo.transparent && nTries < maxTries;
    eventIsOK = !eventInfo.transparent;
    if(eventIsOK) {
      
      for(G4int i = 0; i < eventInfo.nParticles; i++) {
	G4int A = eventInfo.A[i];
	G4int Z = eventInfo.Z[i];
	//	G4cout <<"INCL particle A = " << A << " Z = " << Z << G4endl;
	G4double kinE = eventInfo.EKin[i];
	G4double px = eventInfo.px[i];
	G4double py = eventInfo.py[i];
	G4double pz = eventInfo.pz[i];
	G4DynamicParticle *p = G4INCLXXFactory::toG4Particle(A, Z , kinE, px, py, pz);
	if(p != 0) {
	  const G4LorentzVector momentum = p->Get4Momentum();
	  // Set the four-momentum of the reaction products and apply the toLabFrame rotation 
	  p->Set4Momentum(toLabFrame * momentum);
	  theResult.AddSecondary(p);

	  if(storeDebugOutput) {
	    (*debugOutputFile) << "p " << eventInfo.A[i] << '\t' << eventInfo.Z[i] <<
	      '\t' << eventInfo.emissionTime[i] << '\t' << eventInfo.EKin[i] << '\t' <<
	      eventInfo.px[i] << '\t' << eventInfo.py[i] << '\t' <<
	      eventInfo.pz[i] << '\t' << eventInfo.theta[i] << '\t' <<
	      eventInfo.phi[i] << '\t' << eventInfo.origin[i] << '\t' <<
	      eventInfo.history[i] << std::endl;
	  }
	} else {
	  G4cout <<"Warning: INCL++ produced a particle that couldn't be converted to Geant4 particle." << G4endl;
	}
      }

      for(G4int i = 0; i < eventInfo.nRemnants; i++) {
	G4int A = eventInfo.ARem[i];
	G4int Z = eventInfo.ZRem[i];
	//	G4cout <<"INCL particle A = " << A << " Z = " << Z << G4endl;
	G4double kinE = eventInfo.EKinRem[i];
	G4double px = eventInfo.pxRem[i];
	G4double py = eventInfo.pyRem[i];
	G4double pz = eventInfo.pzRem[i];
	G4double excitationE = eventInfo.EStarRem[i];
	G4double nuclearMass = G4NucleiProperties::GetNuclearMass(A, Z) + excitationE;
	G4double scaling = G4INCLXXFactory::remnant4MomentumScaling(nuclearMass,
								    kinE,
								    px, py, pz);
	G4LorentzVector fourMomentum(scaling * px, scaling * py, scaling * pz,
				     nuclearMass + kinE);
	if(std::abs(scaling - 1.0) > 0.01) {
	  G4cout <<"WARNING: momentum scaling = " << scaling << G4endl;
	  G4cout <<"Lorentz vector = " << fourMomentum << G4endl;
	}
	G4Fragment remnant(A, Z, fourMomentum);
	remnants.push_back(remnant);
	if(storeDebugOutput) {
	  (*debugOutputFile) << "r " << eventInfo.ARem[i] << '\t' <<
	    eventInfo.ZRem[i] << '\t' << eventInfo.EStarRem[i] << '\t' <<
	    eventInfo.JRem[i] << '\t' << eventInfo.EKinRem[i] << '\t' <<
	    eventInfo.pxRem[i] << '\t' << eventInfo.pyRem[i] << '\t' <<
	    eventInfo.pzRem[i] << '\t' << eventInfo.thetaRem[i] << '\t' <<
	    eventInfo.phiRem[i] << std::endl;
	}
      }
    }
    delete theINCLModel;
    nTries++;
  } while(!eventIsOK && nTries < maxTries);
  
  // De-excitation:

  if(theExcitationHandler != 0) {
    for(std::list<G4Fragment>::const_iterator i = remnants.begin();
	i != remnants.end(); i++) {
      const G4LorentzVector remnant4Momentum = (*i).GetMomentum();
      G4LorentzRotation toRemnantZ;
      toRemnantZ.rotateZ(-remnant4Momentum.theta());
      toRemnantZ.rotateY(-remnant4Momentum.phi());
      const G4LorentzRotation toRemnantLab = toRemnantZ.inverse();

      G4LorentzVector remnant4MomentumCM = remnant4Momentum;
      remnant4MomentumCM *= toRemnantZ;
      remnant4MomentumCM.boost(-remnant4Momentum.boostVector());

      G4ReactionProductVector *deExcitationResult = theExcitationHandler->BreakItUp((*i));
    
      for(G4ReactionProductVector::iterator fragment = deExcitationResult->begin();
	  fragment != deExcitationResult->end(); ++fragment) {
	G4ParticleDefinition *def = (*fragment)->GetDefinition();
	if(def != 0) {
	  G4DynamicParticle *theFragment = new G4DynamicParticle(def, (*fragment)->GetMomentum());
	  G4LorentzVector labMomentum = theFragment->Get4Momentum();
	  labMomentum.boost(remnant4Momentum.boostVector());
	  labMomentum *= toRemnantLab;
	  labMomentum *= toLabFrame;
	  theFragment->Set4Momentum(labMomentum);
	  theResult.AddSecondary(theFragment);

	  if(storeDebugOutput) {
	    G4int A = def->GetAtomicMass();
	    G4int Z = def->GetAtomicNumber();
	    G4double fragmentEkin = theFragment->GetKineticEnergy() / MeV;
	    G4ThreeVector mom = theFragment->GetMomentum();
	    (*debugOutputFile) << "de-excitation: p " << A << '\t' << Z << '\t' <<
	      "-1.0" << '\t' << fragmentEkin << '\t' <<
	      mom.x() << '\t' << mom.y() << '\t' << mom.z() << '\t' <<
	      mom.theta() << '\t' <<
	      mom.phi() << '\t' << "-2" << '\t' <<
	      "0" << std::endl;
	  }
	}
      }

      for(G4ReactionProductVector::iterator fragment = deExcitationResult->begin();
	  fragment != deExcitationResult->end(); ++fragment) {
	delete (*fragment);
      }
      deExcitationResult->clear();
      delete deExcitationResult;    
    }
  }

  remnants.clear();

  return &theResult;
}
  
G4ReactionProductVector* G4INCLXXInterface::Propagate(G4KineticTrackVector* , G4V3DNucleus* ) {
  return 0;
}


