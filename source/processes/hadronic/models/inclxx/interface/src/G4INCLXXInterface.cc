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
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include <cmath>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4INCLXXInterface.hh"
#include "G4GenericIon.hh"
#include "G4INCLCascade.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4HadSecondary.hh"
#include "G4ParticleTable.hh"
#include "G4INCLXXInterfaceStore.hh"
#include "G4INCLXXVInterfaceTally.hh"
#include "G4String.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4INCLVersion.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4CompetitiveFission.hh"
#include "G4FissionLevelDensityParameterINCLXX.hh"
#include "G4PhysicsModelCatalog.hh"

#include "G4HyperNucleiProperties.hh"
#include "G4HyperTriton.hh"
#include "G4HyperH4.hh"
#include "G4HyperAlpha.hh"
#include "G4DoubleHyperH4.hh"
#include "G4DoubleHyperDoubleNeutron.hh"
#include "G4HyperHe5.hh"

G4INCLXXInterface::G4INCLXXInterface(G4VPreCompoundModel * const aPreCompound) :
  G4VIntraNuclearTransportModel(G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName()),
  theINCLModel(NULL),
  thePreCompoundModel(aPreCompound),
  theInterfaceStore(G4INCLXXInterfaceStore::GetInstance()),
  theTally(NULL),
  complainedAboutBackupModel(false),
  complainedAboutPreCompound(false),
  theIonTable(G4IonTable::GetIonTable()),
  theINCLXXLevelDensity(NULL),
  theINCLXXFissionProbability(NULL),
  secID(-1)
{
  if(!thePreCompoundModel) {
    G4HadronicInteraction* p =
      G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
    thePreCompoundModel = static_cast<G4VPreCompoundModel*>(p);
    if(!thePreCompoundModel) { thePreCompoundModel = new G4PreCompoundModel; }
  }

  // Use the environment variable G4INCLXX_NO_DE_EXCITATION to disable de-excitation
  if(std::getenv("G4INCLXX_NO_DE_EXCITATION")) {
    G4String message = "de-excitation is completely disabled!";
    theInterfaceStore->EmitWarning(message);
    theDeExcitation = 0;
  } else {
    G4HadronicInteraction* p =
      G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
    theDeExcitation = static_cast<G4VPreCompoundModel*>(p);
    if(!theDeExcitation) { theDeExcitation = new G4PreCompoundModel; }

    // set the fission parameters for G4ExcitationHandler
    G4VEvaporationChannel * const theFissionChannel =
      theDeExcitation->GetExcitationHandler()->GetEvaporation()->GetFissionChannel();
    G4CompetitiveFission * const theFissionChannelCast = dynamic_cast<G4CompetitiveFission *>(theFissionChannel);
    if(theFissionChannelCast) {
      theINCLXXLevelDensity = new G4FissionLevelDensityParameterINCLXX;
      theFissionChannelCast->SetLevelDensityParameter(theINCLXXLevelDensity);
      theINCLXXFissionProbability = new G4FissionProbability;
      theINCLXXFissionProbability->SetFissionLevelDensityParameter(theINCLXXLevelDensity);
      theFissionChannelCast->SetEmissionStrategy(theINCLXXFissionProbability);
      theInterfaceStore->EmitBigWarning("INCL++/G4ExcitationHandler uses its own level-density parameter for fission");
    } else {
      theInterfaceStore->EmitBigWarning("INCL++/G4ExcitationHandler could not use its own level-density parameter for fission");
    }
  }

  // use the envvar G4INCLXX_DUMP_REMNANT to dump information about the
  // remnants on stdout
  if(std::getenv("G4INCLXX_DUMP_REMNANT"))
    dumpRemnantInfo = true;
  else
    dumpRemnantInfo = false;

  theBackupModel = new G4BinaryLightIonReaction;
  theBackupModelNucleon = new G4BinaryCascade;
  secID = G4PhysicsModelCatalog::GetModelID( "model_INCLXXCascade" );
}

G4INCLXXInterface::~G4INCLXXInterface()
{
  delete theINCLXXLevelDensity;
  delete theINCLXXFissionProbability;
}

G4bool G4INCLXXInterface::AccurateProjectile(const G4HadProjectile &aTrack, const G4Nucleus &theNucleus) const {
  // Use direct kinematics if the projectile is a nucleon or a pion
  const G4ParticleDefinition *projectileDef = aTrack.GetDefinition();
  if(std::abs(projectileDef->GetBaryonNumber()) < 2) // Every non-composite particle. abs()-> in case of anti-nucleus projectile
    return false;

  // Here all projectiles should be nuclei
  const G4int pA = projectileDef->GetAtomicMass();
  if(pA<=0) {
    std::stringstream ss;
    ss << "the model does not know how to handle a collision between a "
      << projectileDef->GetParticleName() << " projectile and a Z="
      << theNucleus.GetZ_asInt() << ", A=" << theNucleus.GetA_asInt();
    theInterfaceStore->EmitBigWarning(ss.str());
    return true;
  }

  // If either nucleus is a LCP (A<=4), run the collision as light on heavy
  const G4int tA = theNucleus.GetA_asInt();
  if(tA<=4 || pA<=4) {
    if(pA<tA)
      return false;
    else
      return true;
  }

  // If one of the nuclei is heavier than theMaxProjMassINCL, run the collision
  // as light on heavy.
  // Note that here we are sure that either the projectile or the target is
  // smaller than theMaxProjMass; otherwise theBackupModel would have been
  // called.
  const G4int theMaxProjMassINCL = theInterfaceStore->GetMaxProjMassINCL();
  if(pA > theMaxProjMassINCL)
    return true;
  else if(tA > theMaxProjMassINCL)
    return false;
  else
    // In all other cases, use the global setting
    return theInterfaceStore->GetAccurateProjectile();
}

G4HadFinalState* G4INCLXXInterface::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  G4ParticleDefinition const * const trackDefinition = aTrack.GetDefinition();
  const G4bool isIonTrack = trackDefinition->GetParticleType()==G4GenericIon::GenericIon()->GetParticleType();
  const G4int trackA = trackDefinition->GetAtomicMass();
  const G4int trackZ = (G4int) trackDefinition->GetPDGCharge();
  const G4int trackL = trackDefinition->GetNumberOfLambdasInHypernucleus();
  const G4int nucleusA = theNucleus.GetA_asInt();
  const G4int nucleusZ = theNucleus.GetZ_asInt();

  // For reactions induced by weird projectiles (e.g. He2), bail out
  if((isIonTrack && ((trackZ<=0 && trackL==0) || trackA<=trackZ)) ||
                    (nucleusA>1 && (nucleusZ<=0 || nucleusA<=nucleusZ))) {
    theResult.Clear();
    theResult.SetStatusChange(isAlive);
    theResult.SetEnergyChange(aTrack.GetKineticEnergy());
    theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theResult;
  }

  // For reactions on nucleons, use the backup model (without complaining)
  if(trackA<=1 && nucleusA<=1) {
    return theBackupModelNucleon->ApplyYourself(aTrack, theNucleus);
  }

  // For systems heavier than theMaxProjMassINCL, use another model (typically
  // BIC)
  const G4int theMaxProjMassINCL = theInterfaceStore->GetMaxProjMassINCL();
  if(trackA > theMaxProjMassINCL && nucleusA > theMaxProjMassINCL) {
    if(!complainedAboutBackupModel) {
      complainedAboutBackupModel = true;
      std::stringstream ss;
      ss << "INCL++ refuses to handle reactions between nuclei with A>"
        << theMaxProjMassINCL
        << ". A backup model ("
        << theBackupModel->GetModelName()
        << ") will be used instead.";
      theInterfaceStore->EmitBigWarning(ss.str());
    }
    return theBackupModel->ApplyYourself(aTrack, theNucleus);
  }

  // For energies lower than cascadeMinEnergyPerNucleon, use PreCompound
  const G4double cascadeMinEnergyPerNucleon = theInterfaceStore->GetCascadeMinEnergyPerNucleon();
  const G4double trackKinE = aTrack.GetKineticEnergy();
  if((trackDefinition==G4Neutron::NeutronDefinition() || trackDefinition==G4Proton::ProtonDefinition())
      && trackKinE < cascadeMinEnergyPerNucleon) {
    if(!complainedAboutPreCompound) {
      complainedAboutPreCompound = true;
      std::stringstream ss;
      ss << "INCL++ refuses to handle nucleon-induced reactions below "
        << cascadeMinEnergyPerNucleon / MeV
        << " MeV. A PreCoumpound model ("
        << thePreCompoundModel->GetModelName()
        << ") will be used instead.";
      theInterfaceStore->EmitBigWarning(ss.str());
    }
    return thePreCompoundModel->ApplyYourself(aTrack, theNucleus);
  }

  // Calculate the total four-momentum in the entrance channel
  const G4double theNucleusMass = theIonTable->GetIonMass(nucleusZ, nucleusA);
  const G4double theTrackMass = trackDefinition->GetPDGMass();
  const G4double theTrackEnergy = trackKinE + theTrackMass;
  const G4double theTrackMomentumAbs2 = theTrackEnergy*theTrackEnergy - theTrackMass*theTrackMass;
  const G4double theTrackMomentumAbs = ((theTrackMomentumAbs2>0.0) ? std::sqrt(theTrackMomentumAbs2) : 0.0);
  const G4ThreeVector theTrackMomentum = aTrack.Get4Momentum().getV().unit() * theTrackMomentumAbs;

  G4LorentzVector goodTrack4Momentum(theTrackMomentum, theTrackEnergy);
  G4LorentzVector fourMomentumIn;
  fourMomentumIn.setE(theTrackEnergy + theNucleusMass);
  fourMomentumIn.setVect(theTrackMomentum);

  // Check if inverse kinematics should be used
  const G4bool inverseKinematics = AccurateProjectile(aTrack, theNucleus);

  // If we are running in inverse kinematics, redefine aTrack and theNucleus
  G4LorentzRotation *toInverseKinematics = NULL;
  G4LorentzRotation *toDirectKinematics = NULL;
  G4HadProjectile const *aProjectileTrack = &aTrack;
  G4Nucleus *theTargetNucleus = &theNucleus;
  if(inverseKinematics) {
    G4ParticleDefinition *oldTargetDef = theIonTable->GetIon(nucleusZ, nucleusA, 0);
    const G4ParticleDefinition *oldProjectileDef = trackDefinition;

    if(oldProjectileDef != 0 && oldTargetDef != 0) {
      const G4int newTargetA = oldProjectileDef->GetAtomicMass();
      const G4int newTargetZ = oldProjectileDef->GetAtomicNumber();
      const G4int newTargetL = oldProjectileDef->GetNumberOfLambdasInHypernucleus();
      if(newTargetA > 0 && newTargetZ > 0) {
        // This should give us the same energy per nucleon
        theTargetNucleus = new G4Nucleus(newTargetA, newTargetZ, newTargetL);
        toInverseKinematics = new G4LorentzRotation(goodTrack4Momentum.boostVector());
        G4LorentzVector theProjectile4Momentum(0.0, 0.0, 0.0, theNucleusMass);
        G4DynamicParticle swappedProjectileParticle(oldTargetDef, (*toInverseKinematics) * theProjectile4Momentum);
        aProjectileTrack = new G4HadProjectile(swappedProjectileParticle);
      } else {
        G4String message = "badly defined target after swapping. Falling back to normal (non-swapped) mode.";
        theInterfaceStore->EmitWarning(message);
        toInverseKinematics = new G4LorentzRotation;
      }
    } else {
      G4String message = "oldProjectileDef or oldTargetDef was null";
      theInterfaceStore->EmitWarning(message);
      toInverseKinematics = new G4LorentzRotation;
    }
  }

  // INCL assumes the projectile particle is going in the direction of
  // the Z-axis. Here we construct proper rotation to convert the
  // momentum vectors of the outcoming particles to the original
  // coordinate system.
  G4LorentzVector projectileMomentum = aProjectileTrack->Get4Momentum();

  // INCL++ assumes that the projectile is going in the direction of
  // the z-axis. In principle, if the coordinate system used by G4
  // hadronic framework is defined differently we need a rotation to
  // transform the INCL++ reaction products to the appropriate
  // frame. Please note that it isn't necessary to apply this
  // transform to the projectile because when creating the INCL++
  // projectile particle; \see{toINCLKineticEnergy} needs to use only the
  // projectile energy (direction is simply assumed to be along z-axis).
  G4RotationMatrix toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4RotationMatrix toLabFrame3 = toZ.inverse();
  G4LorentzRotation toLabFrame(toLabFrame3);
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

  const G4int maxTries = 200;
  G4int nTries = 0;
  // INCL can generate transparent events. However, this is meaningful
  // only in the standalone code. In Geant4 we must "force" INCL to
  // produce a valid cascade.
  G4bool eventIsOK = false;
  do {
    const G4INCL::ParticleSpecies theSpecies = toINCLParticleSpecies(*aProjectileTrack);
    const G4double kineticEnergy = toINCLKineticEnergy(*aProjectileTrack);

    // The INCL model will be created at the first use
    theINCLModel = G4INCLXXInterfaceStore::GetInstance()->GetINCLModel();

    const G4INCL::EventInfo eventInfo = theINCLModel->processEvent(theSpecies, kineticEnergy,
								   theTargetNucleus->GetA_asInt(),
								   theTargetNucleus->GetZ_asInt(),
								   -theTargetNucleus->GetL());  // Strangeness has opposite sign
    //    eventIsOK = !eventInfo.transparent && nTries < maxTries;                              // of the number of Lambdas 
    eventIsOK = !eventInfo.transparent;
    if(eventIsOK) {

      // If the collision was run in inverse kinematics, prepare to boost back
      // all the reaction products
      if(inverseKinematics) {
        toDirectKinematics = new G4LorentzRotation(toInverseKinematics->inverse());
      }

      G4LorentzVector fourMomentumOut;

      for(G4int i = 0; i < eventInfo.nParticles; ++i) {
	G4int A = eventInfo.A[i];
        G4int Z = eventInfo.Z[i];
	G4int S = eventInfo.S[i];  // Strangeness
        G4int PDGCode = eventInfo.PDGCode[i];
	//	G4cout <<"INCL particle A = " << A << " Z = " << Z << " S = " << S << G4endl;
	G4double kinE = eventInfo.EKin[i];
	G4double px = eventInfo.px[i];
	G4double py = eventInfo.py[i];
	G4double pz = eventInfo.pz[i];
	G4DynamicParticle *p = toG4Particle(A, Z, S, PDGCode, kinE, px, py, pz);
	if(p != 0) {
	  G4LorentzVector momentum = p->Get4Momentum();

          // Apply the toLabFrame rotation
          momentum *= toLabFrame;
          // Apply the inverse-kinematics boost
          if(inverseKinematics) {
            momentum *= *toDirectKinematics;
            momentum.setVect(-momentum.vect());
          }

	  // Set the four-momentum of the reaction products
	  p->Set4Momentum(momentum);
          fourMomentumOut += momentum;

	  // Propagate the particle's parent resonance information
	  G4HadSecondary secondary(p, 1.0, secID);
	  G4ParticleDefinition* parentResonanceDef = nullptr;
	  if ( eventInfo.parentResonancePDGCode[i] != 0 ) {
	    parentResonanceDef = G4ParticleTable::GetParticleTable()->FindParticle(eventInfo.parentResonancePDGCode[i]);
	  }
	  secondary.SetParentResonanceDef(parentResonanceDef);
	  secondary.SetParentResonanceID(eventInfo.parentResonanceID[i]);
	  
	  theResult.AddSecondary(secondary);

	} else {
	  G4String message = "the model produced a particle that couldn't be converted to Geant4 particle.";
          theInterfaceStore->EmitWarning(message);
	}
      }

      for(G4int i = 0; i < eventInfo.nRemnants; ++i) {
	const G4int A = eventInfo.ARem[i];
	const G4int Z = eventInfo.ZRem[i];
	const G4int S = eventInfo.SRem[i];
	//	G4cout <<"INCL particle A = " << A << " Z = " << Z << " S= " << S << G4endl;
        // Check that the remnant is a physical bound state: if not, resample the collision.
        if(( Z == 0  &&  S == 0  &&  A > 1 ) ||                 // No bound states for nn, nnn, nnnn, ...
           ( Z == 0  &&  S != 0  &&  A < 4 ) ||                 // No bound states for nl, ll, nnl, nll, lll
           ( Z != 0  &&  S != 0  &&  A == Z + std::abs(S) )) {  // No bound states for pl, ppl, pll, ...
	  std::stringstream ss;
	  ss << "unphysical residual fragment : Z=" << Z << "  S=" << S << "  A=" << A 
             << "  skipping it and resampling the collision";
	  theInterfaceStore->EmitWarning(ss.str());
	  eventIsOK = false;
          continue;
	}
	const G4double kinE = eventInfo.EKinRem[i];
	const G4double px = eventInfo.pxRem[i];
	const G4double py = eventInfo.pyRem[i];
	const G4double pz = eventInfo.pzRem[i];
        G4ThreeVector spin(
            eventInfo.jxRem[i]*hbar_Planck,
            eventInfo.jyRem[i]*hbar_Planck,
            eventInfo.jzRem[i]*hbar_Planck
            );
	const G4double excitationE = eventInfo.EStarRem[i];
	G4double nuclearMass = excitationE;
	if ( S == 0 ) {
	  nuclearMass += G4NucleiProperties::GetNuclearMass(A, Z);
	} else {
	  // Assumed that the opposite of the strangeness of the remnant gives the number of Lambdas inside it
	  nuclearMass += G4HyperNucleiProperties::GetNuclearMass(A, Z, std::abs(S));
	}
        const G4double scaling = remnant4MomentumScaling(nuclearMass, kinE, px, py, pz);
	G4LorentzVector fourMomentum(scaling * px, scaling * py, scaling * pz,
				     nuclearMass + kinE);
	if(std::abs(scaling - 1.0) > 0.01) {
          std::stringstream ss;
          ss << "momentum scaling = " << scaling
             << "\n                Lorentz vector = " << fourMomentum
	     << ")\n                A = " << A << ", Z = " << Z << ", S = " << S
             << "\n                E* = " << excitationE << ", nuclearMass = " << nuclearMass
             << "\n                remnant i=" << i << ", nRemnants=" << eventInfo.nRemnants
             << "\n                Reaction was: " << aTrack.GetKineticEnergy()/MeV
             << "-MeV " << trackDefinition->GetParticleName() << " + "
             << theIonTable->GetIonName(theNucleus.GetZ_asInt(), theNucleus.GetA_asInt(), 0)
             << ", in " << (inverseKinematics ? "inverse" : "direct") << " kinematics.";
          theInterfaceStore->EmitWarning(ss.str());
	}

        // Apply the toLabFrame rotation
        fourMomentum *= toLabFrame;
        spin *= toLabFrame3;
        // Apply the inverse-kinematics boost
        if(inverseKinematics) {
          fourMomentum *= *toDirectKinematics;
          fourMomentum.setVect(-fourMomentum.vect());
        }

        fourMomentumOut += fourMomentum;
	G4Fragment remnant(A, Z, std::abs(S), fourMomentum);  // Assumed that -strangeness gives the number of Lambdas
        remnant.SetAngularMomentum(spin);
	remnant.SetCreatorModelID(secID);
        if(dumpRemnantInfo) {
          G4cerr << "G4INCLXX_DUMP_REMNANT: " << remnant << "  spin: " << spin << G4endl;
        }
	remnants.push_back(remnant);
      }

      // Give up is the event is not ok (e.g. unphysical residual)
      if(!eventIsOK) {
        const G4int nSecondaries = (G4int)theResult.GetNumberOfSecondaries();
        for(G4int j=0; j<nSecondaries; ++j) delete theResult.GetSecondary(j)->GetParticle();
        theResult.Clear();
        theResult.SetStatusChange(stopAndKill);
        remnants.clear();
      } else {
        // Check four-momentum conservation
        const G4LorentzVector violation4Momentum = fourMomentumOut - fourMomentumIn;
        const G4double energyViolation = std::abs(violation4Momentum.e());
        const G4double momentumViolation = violation4Momentum.rho();
        if(energyViolation > G4INCLXXInterfaceStore::GetInstance()->GetConservationTolerance()) {
          std::stringstream ss;
          ss << "energy conservation violated by " << energyViolation/MeV << " MeV in "
             << aTrack.GetKineticEnergy()/MeV << "-MeV " << trackDefinition->GetParticleName()
             << " + " << theIonTable->GetIonName(theNucleus.GetZ_asInt(), theNucleus.GetA_asInt(), 0)
             << " inelastic reaction, in " << (inverseKinematics ? "inverse" : "direct") << " kinematics. Will resample.";
          theInterfaceStore->EmitWarning(ss.str());
          eventIsOK = false;
          const G4int nSecondaries = (G4int)theResult.GetNumberOfSecondaries();
          for(G4int j=0; j<nSecondaries; ++j) delete theResult.GetSecondary(j)->GetParticle();
          theResult.Clear();
          theResult.SetStatusChange(stopAndKill);
          remnants.clear();
        } else if(momentumViolation > G4INCLXXInterfaceStore::GetInstance()->GetConservationTolerance()) {
          std::stringstream ss;
          ss << "momentum conservation violated by " << momentumViolation/MeV << " MeV in "
             << aTrack.GetKineticEnergy()/MeV << "-MeV " << trackDefinition->GetParticleName()
             << " + " << theIonTable->GetIonName(theNucleus.GetZ_asInt(), theNucleus.GetA_asInt(), 0)
             << " inelastic reaction, in " << (inverseKinematics ? "inverse" : "direct") << " kinematics. Will resample.";
          theInterfaceStore->EmitWarning(ss.str());
          eventIsOK = false;
          const G4int nSecondaries = (G4int)theResult.GetNumberOfSecondaries();
          for(G4int j=0; j<nSecondaries; ++j) delete theResult.GetSecondary(j)->GetParticle();
          theResult.Clear();
          theResult.SetStatusChange(stopAndKill);
          remnants.clear();
        }
      }
    }
    nTries++;
  } while(!eventIsOK && nTries < maxTries); /* Loop checking, 10.07.2015, D.Mancusi */

  // Clean up the objects that we created for the inverse kinematics
  if(inverseKinematics) {
    delete aProjectileTrack;
    delete theTargetNucleus;
    delete toInverseKinematics;
    delete toDirectKinematics;
  }

  if(!eventIsOK) {
    std::stringstream ss;
    ss << "maximum number of tries exceeded for the proposed "
    << aTrack.GetKineticEnergy()/MeV << "-MeV " << trackDefinition->GetParticleName()
    << " + " << theIonTable->GetIonName(nucleusZ, nucleusA, 0)
    << " inelastic reaction, in " << (inverseKinematics ? "inverse" : "direct") << " kinematics.";
    theInterfaceStore->EmitWarning(ss.str());
    theResult.SetStatusChange(isAlive);
    theResult.SetEnergyChange(aTrack.GetKineticEnergy());
    theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theResult;
  }

  // De-excitation:

  if(theDeExcitation != 0) {
    for(std::list<G4Fragment>::iterator i = remnants.begin();
	i != remnants.end(); ++i) {
      G4ReactionProductVector *deExcitationResult = theDeExcitation->DeExcite((*i));

      for(G4ReactionProductVector::iterator fragment = deExcitationResult->begin();
	  fragment != deExcitationResult->end(); ++fragment) {
	const G4ParticleDefinition *def = (*fragment)->GetDefinition();
	if(def != 0) {
	  G4DynamicParticle *theFragment = new G4DynamicParticle(def, (*fragment)->GetMomentum());
	  theResult.AddSecondary(theFragment, (*fragment)->GetCreatorModelID());
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

  if((theTally = theInterfaceStore->GetTally()))
    theTally->Tally(aTrack, theNucleus, theResult);

  return &theResult;
}

G4ReactionProductVector* G4INCLXXInterface::Propagate(G4KineticTrackVector* , G4V3DNucleus* ) {
  return 0;
}

G4INCL::ParticleType G4INCLXXInterface::toINCLParticleType(G4ParticleDefinition const * const pdef) const {
  if(     pdef == G4Proton::Proton())               return G4INCL::Proton;
  else if(pdef == G4Neutron::Neutron())             return G4INCL::Neutron;
  else if(pdef == G4PionPlus::PionPlus())           return G4INCL::PiPlus;
  else if(pdef == G4PionMinus::PionMinus())         return G4INCL::PiMinus;
  else if(pdef == G4PionZero::PionZero())           return G4INCL::PiZero;
  else if(pdef == G4KaonPlus::KaonPlus())           return G4INCL::KPlus;
  else if(pdef == G4KaonZero::KaonZero())           return G4INCL::KZero;
  else if(pdef == G4KaonMinus::KaonMinus())         return G4INCL::KMinus;
  else if(pdef == G4AntiKaonZero::AntiKaonZero())   return G4INCL::KZeroBar;
  // For K0L & K0S we do not take into account K0/K0B oscillations
  else if(pdef == G4KaonZeroLong::KaonZeroLong())   return G4UniformRand() < 0.5 ? G4INCL::KZeroBar : G4INCL::KZero;
  else if(pdef == G4KaonZeroShort::KaonZeroShort()) return G4UniformRand() < 0.5 ? G4INCL::KZeroBar : G4INCL::KZero; 
  else if(pdef == G4Deuteron::Deuteron())           return G4INCL::Composite;
  else if(pdef == G4Triton::Triton())               return G4INCL::Composite;
  else if(pdef == G4He3::He3())                     return G4INCL::Composite;
  else if(pdef == G4Alpha::Alpha())                 return G4INCL::Composite;
  else if(pdef->GetParticleType() == G4GenericIon::GenericIon()->GetParticleType()) return G4INCL::Composite;
  else                                              return G4INCL::UnknownParticle;
}

G4INCL::ParticleSpecies G4INCLXXInterface::toINCLParticleSpecies(G4HadProjectile const &aTrack) const {
  const G4ParticleDefinition *pdef = aTrack.GetDefinition();
  const G4INCL::ParticleType theType = toINCLParticleType(pdef);
  if(theType!=G4INCL::Composite)
    return G4INCL::ParticleSpecies(theType);
  else {
    G4INCL::ParticleSpecies theSpecies;
    theSpecies.theType=theType;
    theSpecies.theA=pdef->GetAtomicMass();
    theSpecies.theZ=pdef->GetAtomicNumber();
    return theSpecies;
  }
}

G4double G4INCLXXInterface::toINCLKineticEnergy(G4HadProjectile const &aTrack) const {
  return aTrack.GetKineticEnergy();
}

G4ParticleDefinition *G4INCLXXInterface::toG4ParticleDefinition(G4int A, G4int Z, G4int S, G4int PDGCode) const {
  if       (PDGCode == 2212) { return G4Proton::Proton();
  } else if(PDGCode == 2112) { return G4Neutron::Neutron();
  } else if(PDGCode == 211)  { return G4PionPlus::PionPlus();
  } else if(PDGCode == 111)  { return G4PionZero::PionZero();
  } else if(PDGCode == -211) { return G4PionMinus::PionMinus();
  
  } else if(PDGCode == 221)  { return G4Eta::Eta();
  } else if(PDGCode == 22)   { return G4Gamma::Gamma();
  
  } else if(PDGCode == 3122) { return G4Lambda::Lambda();
  } else if(PDGCode == 3222) { return G4SigmaPlus::SigmaPlus();
  } else if(PDGCode == 3212) { return G4SigmaZero::SigmaZero();
  } else if(PDGCode == 3112) { return G4SigmaMinus::SigmaMinus();
  } else if(PDGCode == 321)  { return G4KaonPlus::KaonPlus();
  } else if(PDGCode == -321) { return G4KaonMinus::KaonMinus();
  } else if(PDGCode == 130)  { return G4KaonZeroLong::KaonZeroLong();
  } else if(PDGCode == 310)  { return G4KaonZeroShort::KaonZeroShort();
  
  } else if(PDGCode == 1002) { return G4Deuteron::Deuteron();
  } else if(PDGCode == 1003) { return G4Triton::Triton();
  } else if(PDGCode == 2003) { return G4He3::He3();
  } else if(PDGCode == 2004) { return G4Alpha::Alpha();
  } else if(S != 0) {  // Assumed that -S gives the number of Lambdas
    if (A == 3 && Z == 1 && S == -1 ) return G4HyperTriton::Definition();
    if (A == 4 && Z == 1 && S == -1 ) return G4HyperH4::Definition();
    if (A == 4 && Z == 2 && S == -1 ) return G4HyperAlpha::Definition();
    if (A == 4 && Z == 1 && S == -2 ) return G4DoubleHyperH4::Definition();
    if (A == 4 && Z == 0 && S == -2 ) return G4DoubleHyperDoubleNeutron::Definition();
    if (A == 5 && Z == 2 && S == -1 ) return G4HyperHe5::Definition();
  } else if(A > 0 && Z > 0 && A > Z) { // Returns ground state ion definition.
    return theIonTable->GetIon(Z, A, 0);
  }
  return 0; // Error, unrecognized particle
}

G4DynamicParticle *G4INCLXXInterface::toG4Particle(G4int A, G4int Z, G4int S, G4int PDGCode,
						   G4double kinE, G4double px,
                                                   G4double py, G4double pz) const {
  const G4ParticleDefinition *def = toG4ParticleDefinition(A, Z, S, PDGCode);
  if(def == 0) { // Check if we have a valid particle definition
    return 0;
  }
  const G4double energy = kinE * MeV;
  const G4ThreeVector momentum(px, py, pz);
  const G4ThreeVector momentumDirection = momentum.unit();
  G4DynamicParticle *p = new G4DynamicParticle(def, momentumDirection, energy);
  return p;
}

G4double G4INCLXXInterface::remnant4MomentumScaling(G4double mass,
						  G4double kineticE,
						  G4double px, G4double py,
						  G4double pz) const {
  const G4double p2 = px*px + py*py + pz*pz;
  if(p2 > 0.0) {
    const G4double pnew2 = kineticE*kineticE + 2.0*kineticE*mass;
    return std::sqrt(pnew2)/std::sqrt(p2);
  } else {
    return 1.0;
  }
}

void G4INCLXXInterface::ModelDescription(std::ostream& outFile) const {
   outFile
     << "The Liège Intranuclear Cascade (INCL++) is a model for reactions induced\n"
     << "by nucleons, pions and light ion on any nucleus. The reaction is\n"
     << "described as an avalanche of binary nucleon-nucleon collisions, which can\n"
     << "lead to the emission of energetic particles and to the formation of an\n"
     << "excited thermalised nucleus (remnant). The de-excitation of the remnant is\n"
     << "outside the scope of INCL++ and is typically described by another model.\n\n"
     << "INCL++ has been reasonably well tested for nucleon (~50 MeV to ~15 GeV),\n"
     << "pion (idem) and light-ion projectiles (up to A=18, ~10A MeV to 1A GeV).\n"
     << "Most tests involved target nuclei close to the stability valley, with\n"
     << "numbers between 4 and 250.\n\n"
     << "Reference: D. Mancusi et al., Phys. Rev. C90 (2014) 054602\n\n";
}

G4String const &G4INCLXXInterface::GetDeExcitationModelName() const {
  return theDeExcitation->GetModelName();
}
