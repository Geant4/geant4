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
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLCascade.cc
 *
 * INCL Cascade
 */
#include "G4INCLCascade.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLGeant4Random.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobalInfo.hh"

#include "G4INCLPauliBlocking.hh"
#include "G4INCLIPauli.hh"
#include "G4INCLPauliStrict.hh"
#include "G4INCLPauliStandard.hh"
#include "G4INCLPauliStrictStandard.hh"
#include "G4INCLPauliGlobal.hh"
#include "G4INCLCDPP.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"

#include "G4INCLCoulombDistortion.hh"
#include "G4INCLICoulomb.hh"
#include "G4INCLCoulombNone.hh"
#include "G4INCLCoulombNonRelativistic.hh"

#include "G4INCLClustering.hh"
#include "G4INCLClusteringModelIntercomparison.hh"
#include "G4INCLClusteringModelNone.hh"

#include "G4INCLIntersection.hh"

#include "G4INCLCrossSections.hh"

#include <cstring>
#include <cstdlib>
#include <numeric>

namespace G4INCL {

  INCL::INCL(G4INCL::Config const * const config)
    :propagationModel(0), theA(208), theZ(82),
    targetInitSuccess(false),
    maxImpactParameter(0.),
    maxUniverseRadius(0.),
    maxInteractionDistance(0.),
    fixedImpactParameter(0.),
    theConfig(config),
    nucleus(NULL),
    minRemnantSize(4)
  {
    // Set the logger object.
    G4INCL::Logger::setLoggerSlave(new G4INCL::LoggerSlave(theConfig->getLogFileName()));
    G4INCL::Logger::setVerbosityLevel(theConfig->getVerbosity());

    // Set the random number generator algorithm. The system can support
    // multiple different generator algorithms in a completely
    // transparent way.
#ifdef INCLXX_IN_GEANT4_MODE
    G4INCL::Random::setGenerator(new G4INCL::Geant4RandomGenerator());
#else
    G4INCL::Random::setGenerator(new G4INCL::Ranecu(theConfig->getRandomSeeds()));
#endif // INCLXX_IN_GEANT4_MODE

    // Select the Pauli blocking algorithm:
    G4INCL::PauliType pauli = theConfig->getPauliType();
    if(pauli == G4INCL::StrictStatisticalPauli)
      G4INCL::Pauli::setBlocker(new G4INCL::PauliStrictStandard);
    else if(pauli == G4INCL::StatisticalPauli)
      G4INCL::Pauli::setBlocker(new G4INCL::PauliStandard);
    else if(pauli == G4INCL::StrictPauli)
      G4INCL::Pauli::setBlocker(new G4INCL::PauliStrict);
    else if(pauli == G4INCL::GlobalPauli)
      G4INCL::Pauli::setBlocker(new G4INCL::PauliGlobal);
    else if(pauli == G4INCL::NoPauli)
      G4INCL::Pauli::setBlocker(NULL);

    if(theConfig->getCDPP())
      G4INCL::Pauli::setCDPP(new G4INCL::CDPP);
    else
      G4INCL::Pauli::setCDPP(NULL);

    // Select the Coulomb-distortion algorithm:
    G4INCL::CoulombType coulombType = theConfig->getCoulombType();
    if(coulombType == G4INCL::NonRelativisticCoulomb)
      G4INCL::CoulombDistortion::setCoulomb(new G4INCL::CoulombNonRelativistic);
    else // if(coulombType == G4INCL::NoCoulomb)
      G4INCL::CoulombDistortion::setCoulomb(new G4INCL::CoulombNone);

    // Select the clustering algorithm:
    G4INCL::ClusterAlgorithmType clusterAlgorithm = theConfig->getClusterAlgorithm();
    if(clusterAlgorithm == G4INCL::IntercomparisonClusterAlgorithm)
      G4INCL::Clustering::setClusteringModel(new G4INCL::ClusteringModelIntercomparison(theConfig));
    else // if(clusterAlgorithm == G4INCL::NoClusterAlgorithm)
      G4INCL::Clustering::setClusteringModel(new G4INCL::ClusteringModelNone);

    // Initialize the INCL particle table:
    G4INCL::ParticleTable::initialize(theConfig);

    // Propagation model is responsible for finding avatars and
    // transporting the particles. In principle this step is "hidden"
    // behind an abstract interface and the rest of the system does not
    // care how the transportation and avatar finding is done. This
    // should allow us to "easily" experiment with different avatar
    // finding schemes and even to support things like curved
    // trajectories in the future.
    propagationModel = new G4INCL::StandardPropagationModel(theConfig->getLocalEnergyBBType(),theConfig->getLocalEnergyPiType());
    eventAction = new EventAction();
    propagationAction = new PropagationAction();
    avatarAction = new AvatarAction();

    theGlobalInfo.cascadeModel = theConfig->getVersionID().c_str();
    theGlobalInfo.deexcitationModel = "none";

#ifndef INCLXX_IN_GEANT4_MODE
    // Fill in the global information
    theGlobalInfo.At = theConfig->getTargetA();
    theGlobalInfo.Zt = theConfig->getTargetZ();
    const ParticleSpecies theSpecies = theConfig->getProjectileSpecies();
    theGlobalInfo.Ap = theSpecies.theA;
    theGlobalInfo.Zp = theSpecies.theZ;
    theGlobalInfo.Ep = theConfig->getProjectileKineticEnergy();
    // Echo the input parameters to the log file
    INFO(theConfig->echo() << std::endl);
#endif

    fixedImpactParameter = theConfig->getImpactParameter();
  }

  INCL::~INCL() {
    G4INCL::Pauli::deleteBlockers();
    G4INCL::CoulombDistortion::deleteCoulomb();
    G4INCL::Random::deleteGenerator();
    G4INCL::Clustering::deleteClusteringModel();
    G4INCL::Logger::deleteLoggerSlave();
    G4INCL::NuclearDensityFactory::clearCache();
    delete avatarAction;
    delete propagationAction;
    delete eventAction;
    delete propagationModel;
    delete theConfig;
  }

  G4bool INCL::prepareReaction(const ParticleSpecies &projectileSpecies, const G4double kineticEnergy, const G4int A, const G4int Z) {
    if(A < 0 || A > 300 || Z < 1 || Z > 200) {
      ERROR("Unsupported target: A = " << A << " Z = " << Z << std::endl);
      ERROR("Target configuration rejected." << std::endl);
      return false;
    }

    // Initialise the maximum universe radius
    initUniverseRadius(projectileSpecies, kineticEnergy, A, Z);

    // Initialise the nucleus
    theZ = Z;
    if(theConfig->isNaturalTarget())
      theA = ParticleTable::drawRandomNaturalIsotope(Z);
    else
      theA = A;
    initializeTarget(theA, theZ);

    // Set the maximum impact parameter
    maxImpactParameter = CoulombDistortion::maxImpactParameter(projectileSpecies, kineticEnergy, nucleus);
    initMaxInteractionDistance(projectileSpecies, kineticEnergy); // for forced CN events

    // Set the geometric cross section
    theGlobalInfo.geometricCrossSection =
      Math::tenPi*std::pow(maxImpactParameter,2);

    // Set the minimum remnant size
    if(projectileSpecies.theA > 0)
      minRemnantSize = std::min(theA, 4);
    else
      minRemnantSize = std::min(theA-1, 4);

    return true;
  }

  G4bool INCL::initializeTarget(const G4int A, const G4int Z) {
    delete nucleus;

    nucleus = new Nucleus(A, Z, theConfig, maxUniverseRadius);
    nucleus->getStore()->getBook()->reset();
    nucleus->initializeParticles();

    propagationModel->setNucleus(nucleus);
    return true;
  }

  const EventInfo &INCL::processEvent(
      ParticleSpecies const &projectileSpecies,
      const G4double kineticEnergy,
      const G4int targetA,
      const G4int targetZ
      ) {
    // Set the target and the projectile
    targetInitSuccess = prepareReaction(projectileSpecies, kineticEnergy, targetA, targetZ);

    if(!targetInitSuccess) {
      WARN("Target initialisation failed for A=" << targetA << ", Z=" << targetZ << std::endl);
      theEventInfo.transparent=true;
      return theEventInfo;
    }

    const G4bool canRunCascade = preCascade(projectileSpecies, kineticEnergy);
    if(canRunCascade) {
      cascade();
      postCascade();
    }
    return theEventInfo;
  }

  G4bool INCL::preCascade(ParticleSpecies const projectileSpecies, const G4double kineticEnergy) {
    // Reset theEventInfo
    theEventInfo.reset();

    EventInfo::eventNumber++;

    // Increment the global counter for the number of shots
    theGlobalInfo.nShots++;

    // Fill in the event information
    theEventInfo.projectileType = projectileSpecies.theType;
    theEventInfo.Ap = projectileSpecies.theA;
    theEventInfo.Zp = projectileSpecies.theZ;
    theEventInfo.Ep = kineticEnergy;
    theEventInfo.At = nucleus->getA();
    theEventInfo.Zt = nucleus->getZ();

    // Do nothing below the Coulomb barrier
    if(maxImpactParameter<=0.) {
      // Increment the global counter for the number of transparents
      theGlobalInfo.nTransparents++;

      // Fill in the event information
      theEventInfo.transparent = true;

      return false;
    }

    // Randomly draw an impact parameter or use a fixed value, depending on the
    // Config option
    G4double impactParameter, phi;
    if(fixedImpactParameter<0.) {
      impactParameter = maxImpactParameter * std::sqrt(Random::shoot0());
      phi = Random::shoot() * Math::twoPi;
    } else {
      impactParameter = fixedImpactParameter;
      phi = 0.;
    }

    // Fill in the event information
    theEventInfo.impactParameter = impactParameter;

    const G4double effectiveImpactParameter = propagationModel->shoot(projectileSpecies, kineticEnergy, impactParameter, phi);
    if(effectiveImpactParameter < 0.) {
      // Increment the global counter for the number of transparents
      theGlobalInfo.nTransparents++;

      // Fill in the event information
      theEventInfo.transparent = true;

      return false;
    }

    // Fill in the event information
    theEventInfo.transparent = false;
    theEventInfo.effectiveImpactParameter = effectiveImpactParameter;

    return true;
  }

  void INCL::cascade() {
    do {
      // Run book keeping actions that should take place before propagation:
      propagationAction->beforePropagationAction(propagationModel);

      // Get the avatar with the smallest time and propagate particles
      // to that point in time.
      G4INCL::IAvatar *avatar = propagationModel->propagate();

      // Run book keeping actions that should take place after propagation:
      propagationAction->afterPropagationAction(propagationModel, avatar);

      if(avatar == 0) break; // No more avatars in the avatar list.

      // Run book keeping actions that should take place before avatar:
      avatarAction->beforeAvatarAction(avatar, nucleus);

      // Channel is responsible for calculating the outcome of the
      // selected avatar. There are different kinds of channels. The
      // class IChannel is, again, an abstract interface that defines
      // the externally observable behavior of all interaction
      // channels.
      // The handling of the channel is transparent to the API.
      // Final state tells what changed...
      G4INCL::FinalState *finalState = avatar->getFinalState();

      // Run book keeping actions that should take place after avatar:
      avatarAction->afterAvatarAction(avatar, nucleus, finalState);

      // So now we must give this information to the nucleus
      nucleus->applyFinalState(finalState);
      // and now we are ready to process the next avatar!

      delete avatar;
      delete finalState;
    } while(continueCascade());

  }

  void INCL::postCascade() {
    // Fill in the event information
    theEventInfo.stoppingTime = propagationModel->getCurrentTime();

    // Forced CN?
    if(nucleus->getTryCompoundNucleus()) {
      DEBUG("Trying compound nucleus" << std::endl);
      makeCompoundNucleus();
      // Global checks of conservation laws
#ifndef INCLXX_IN_GEANT4_MODE
      if(!theEventInfo.transparent) globalConservationChecks(true);
#endif
      return;
    }

    theEventInfo.transparent = nucleus->isEventTransparent();

    if(theEventInfo.transparent) {
      // Increment the global counter for the number of transparents
      theGlobalInfo.nTransparents++;
      if(nucleus->isForcedTransparent())
        theGlobalInfo.nForcedTransparents++;
      ProjectileRemnant * const projectileRemnant = nucleus->getProjectileRemnant();
      if(projectileRemnant) {
        // Delete the projectile remnant and the particles it contains
        projectileRemnant->deleteParticles();
        nucleus->deleteProjectileRemnant();
        nucleus->getStore()->clearIncoming();
      } else {
        // Clean up the incoming list and force a transparent gracefully
        nucleus->getStore()->deleteIncoming();
      }
    } else {
      // Check if the nucleus contains deltas
      theEventInfo.deltasInside = nucleus->containsDeltas();

      // Take care of any remaining deltas
      theEventInfo.forcedDeltasOutside = nucleus->decayOutgoingDeltas();
      theEventInfo.forcedDeltasInside = nucleus->decayInsideDeltas();

      // Apply Coulomb distortion, if appropriate
      // Note that this will apply Coulomb distortion also on pions emitted by
      // unphysical remnants (see decayInsideDeltas). This is at variance with
      // what INCL4.6 does, but these events are (should be!) so rare that
      // whatever we do doesn't (shouldn't!) make any noticeable difference.
      G4INCL::CoulombDistortion::distortOut(nucleus->getStore()->getOutgoingParticles(), nucleus);

      // If the normal cascade predicted complete fusion, use the tabulated
      // masses to compute the excitation energy, the recoil, etc.
      if(nucleus->getStore()->getOutgoingParticles().size()==0
          && nucleus->getProjectileRemnant()
          && nucleus->getProjectileRemnant()->getParticles().size()==0) {

        DEBUG("Cascade resulted in complete fusion, using realistic fusion kinematics" << std::endl);

        nucleus->useFusionKinematics();
        nucleus->deleteProjectileRemnant();

        if(nucleus->getExcitationEnergy()<0.) {
          // Complete fusion is energetically impossible, return a transparent
          WARN("Complete-fusion kinematics yields negative excitation energy, returning a transparent!" << std::endl);
          theEventInfo.transparent = true;
          return;
        }

      } else { // Normal cascade here

        // Set the excitation energy
        nucleus->setExcitationEnergy(nucleus->computeExcitationEnergy());

        // Make a projectile pre-fragment out of the geometrical and dynamical
        // spectators
        theEventInfo.nUnmergedSpectators = makeProjectileRemnant();

        // Compute recoil momentum, energy and spin of the nucleus
        nucleus->computeRecoilKinematics();

#ifndef INCLXX_IN_GEANT4_MODE
        // Global checks of conservation laws
        globalConservationChecks(false);
#endif

        // Make room for the remnant recoil by rescaling the energies of the
        // outgoing particles.
        if(nucleus->hasRemnant()) rescaleOutgoingForRecoil();

      }

      // Cluster decay
      theEventInfo.clusterDecay = nucleus->decayOutgoingClusters() | nucleus->decayMe();

#ifndef INCLXX_IN_GEANT4_MODE
      // Global checks of conservation laws
      globalConservationChecks(true);
#endif

      // Fill the EventInfo structure
      nucleus->fillEventInfo(&theEventInfo);

      // Check if we have an absorption:
      if(theEventInfo.nucleonAbsorption) theGlobalInfo.nNucleonAbsorptions++;
      if(theEventInfo.pionAbsorption) theGlobalInfo.nPionAbsorptions++;
    }
  }

  void INCL::makeCompoundNucleus() {
    // Reset the internal Nucleus variables
    nucleus->getStore()->clearIncoming();
    nucleus->getStore()->clearOutgoing();
    nucleus->getProjectileRemnant()->reset();
    nucleus->setA(theEventInfo.At);
    nucleus->setZ(theEventInfo.Zt);

    // CN kinematical variables
    // Note: the CN orbital angular momentum is neglected in what follows. We
    // should actually take it into account!
    ThreeVector theCNMomentum = nucleus->getIncomingMomentum();
    ThreeVector theCNSpin = nucleus->getIncomingAngularMomentum();
    const G4double theTargetMass = ParticleTable::getTableMass(theEventInfo.At, theEventInfo.Zt);
    G4int theCNA=theEventInfo.At, theCNZ=theEventInfo.Zt;
    Cluster * const theProjectileRemnant = nucleus->getProjectileRemnant();
    G4double theCNEnergy = theTargetMass + theProjectileRemnant->getEnergy();

    // Loop over the potential participants
    ParticleList initialProjectileComponents = theProjectileRemnant->getParticles();
    std::vector<Particle *> shuffledComponents(initialProjectileComponents.begin(), initialProjectileComponents.end());
    // Shuffle the list of potential participants
    std::random_shuffle(shuffledComponents.begin(), shuffledComponents.end(), shuffleComponentsHelper);

    G4bool success = true;
    G4bool atLeastOneNucleonEntering = false;
    for(std::vector<Particle*>::const_iterator p=shuffledComponents.begin(); p!=shuffledComponents.end(); ++p) {
      // Skip geometrical spectators
      Intersection intersectionUniverse(IntersectionFactory::getEarlierTrajectoryIntersection(
            (*p)->getPosition(),
            (*p)->getPropagationVelocity(),
            maxUniverseRadius));
      if(!intersectionUniverse.exists)
        continue;

      // At least one nucleon must enter the interaction distance
      Intersection intersectionInteraction(IntersectionFactory::getEarlierTrajectoryIntersection(
            (*p)->getPosition(),
            (*p)->getPropagationVelocity(),
            maxInteractionDistance));
      if(intersectionInteraction.exists)
        atLeastOneNucleonEntering = true;

      // Build an entry avatar for this nucleon
      ParticleEntryAvatar theAvatar(0.0, nucleus, *p);
      FinalState *fs = theAvatar.getFinalState();
      nucleus->applyFinalState(fs);
      FinalStateValidity validity = fs->getValidity();
      delete fs;
      switch(validity) {
        case ValidFS:
        case ParticleBelowFermiFS:
          // Add the particle to the CN
          theCNA++;
          theCNZ += (*p)->getZ();
          break;
        case ParticleBelowZeroFS:
        case PauliBlockedFS:
        case NoEnergyConservationFS:
        default:
          success = false;
          break;
      }
    }

    if(!success || !atLeastOneNucleonEntering) {
      DEBUG("No nucleon entering in forced CN, or some nucleons entering below zero, forcing a transparent" << std::endl);
      theEventInfo.transparent = true;
      theGlobalInfo.nTransparents++;
      theGlobalInfo.nForcedTransparents++;
      nucleus->getProjectileRemnant()->deleteParticles();
      nucleus->deleteProjectileRemnant();
      return;
    }

// assert(theCNA==nucleus->getA());
// assert(theCNA>theEventInfo.At);

    // Update the kinematics of the CN
    theCNEnergy -= theProjectileRemnant->getEnergy();
    theCNMomentum -= theProjectileRemnant->getMomentum();

    // Deal with the projectile remnant
    nucleus->finalizeProjectileRemnant(propagationModel->getCurrentTime());

    // Subtract the angular momentum of the projectile remnant
    ParticleList const &outgoing = nucleus->getStore()->getOutgoingParticles();
// assert(outgoing.size()==0 || outgoing.size()==1);
    for(ParticleIter i=outgoing.begin(); i!=outgoing.end(); ++i) {
      theCNSpin -= (*i)->getAngularMomentum();
    }

    // Compute the excitation energy of the CN
    const G4double theCNMass = ParticleTable::getTableMass(theCNA,theCNZ);
    const G4double theCNInvariantMassSquared = theCNEnergy*theCNEnergy-theCNMomentum.mag2();
    if(theCNInvariantMassSquared<0.) {
      // Negative invariant mass squared, return a transparent
      theGlobalInfo.nTransparents++;
      theGlobalInfo.nForcedTransparents++;
      theEventInfo.transparent = true;
      return;
    }
    const G4double theCNExcitationEnergy = std::sqrt(theCNInvariantMassSquared) - theCNMass;
    if(theCNExcitationEnergy<0.) {
      // Negative excitation energy, return a transparent
      theGlobalInfo.nTransparents++;
      theGlobalInfo.nForcedTransparents++;
      theEventInfo.transparent = true;
      return;
    } else {
      // Positive excitation energy, can make a CN
      nucleus->setA(theCNA);
      nucleus->setZ(theCNZ);
      nucleus->setMomentum(theCNMomentum);
      nucleus->setEnergy(theCNEnergy);
      nucleus->setExcitationEnergy(theCNExcitationEnergy);
      nucleus->setMass(theCNMass+theCNExcitationEnergy);
      nucleus->setSpin(theCNSpin); // neglects any orbital angular momentum of the CN

      // Take care of any remaining deltas
      theEventInfo.forcedDeltasOutside = nucleus->decayOutgoingDeltas();

      // Cluster decay
      theEventInfo.clusterDecay = nucleus->decayOutgoingClusters() | nucleus->decayMe();

      // Fill the EventInfo structure
      nucleus->fillEventInfo(&theEventInfo);
      theGlobalInfo.nForcedCompoundNucleus++;
    }
  }

  void INCL::rescaleOutgoingForRecoil() {
    RecoilCMFunctor theRecoilFunctor(nucleus, theEventInfo);

    // Apply the root-finding algorithm
    const G4bool success = RootFinder::solve(&theRecoilFunctor, 1.0);
    if(success) {
      std::pair<G4double,G4double> theSolution = RootFinder::getSolution();
      theRecoilFunctor(theSolution.first); // Apply the solution
    } else {
      WARN("Couldn't accommodate remnant recoil while satisfying energy conservation, root-finding algorithm failed." << std::endl);
    }

  }

#ifndef INCLXX_IN_GEANT4_MODE
  void INCL::globalConservationChecks(G4bool afterRecoil) {
    Nucleus::ConservationBalance theBalance = nucleus->getConservationBalance(theEventInfo,afterRecoil);

    // Global conservation checks
    const G4double pLongBalance = theBalance.momentum.getZ();
    const G4double pTransBalance = theBalance.momentum.perp();
    if(theBalance.Z != 0) {
      ERROR("Violation of charge conservation! ZBalance = " << theBalance.Z << std::endl);
    }
    if(theBalance.A != 0) {
      ERROR("Violation of baryon-number conservation! ABalance = " << theBalance.A << std::endl);
    }
    G4double EThreshold, pLongThreshold, pTransThreshold;
    if(afterRecoil) {
      // Less stringent checks after accommodating recoil
      EThreshold = 10.; // MeV
      pLongThreshold = 1.; // MeV/c
      pTransThreshold = 1.; // MeV/c
    } else {
      // More stringent checks before accommodating recoil
      EThreshold = 0.1; // MeV
      pLongThreshold = 0.1; // MeV/c
      pTransThreshold = 0.1; // MeV/c
    }
    if(std::abs(theBalance.energy)>EThreshold) {
      WARN("Violation of energy conservation > " << EThreshold << " MeV. EBalance = " << theBalance.energy << " afterRecoil = " << afterRecoil << " eventNumber=" << theEventInfo.eventNumber << std::endl);
    }
    if(std::abs(pLongBalance)>pLongThreshold) {
      WARN("Violation of longitudinal momentum conservation > " << pLongThreshold << " MeV/c. pLongBalance = " << pLongBalance << " afterRecoil = " << afterRecoil << " eventNumber=" << theEventInfo.eventNumber << std::endl);
    }
    if(std::abs(pTransBalance)>pTransThreshold) {
      WARN("Violation of transverse momentum conservation > " << pTransThreshold << " MeV/c. pTransBalance = " << pTransBalance << " afterRecoil = " << afterRecoil << " eventNumber=" << theEventInfo.eventNumber << std::endl);
    }

    // Feed the EventInfo variables
    theEventInfo.EBalance = theBalance.energy;
    theEventInfo.pLongBalance = pLongBalance;
    theEventInfo.pTransBalance = pTransBalance;
  }
#endif

  G4bool INCL::continueCascade() {
    // Stop if we have passed the stopping time
    if(propagationModel->getCurrentTime() > propagationModel->getStoppingTime()) {
      DEBUG("Cascade time (" << propagationModel->getCurrentTime()
          << ") exceeded stopping time (" << propagationModel->getStoppingTime()
          << "), stopping cascade" << std::endl);
      return false;
    }
    // Stop if there are no participants and no pions inside the nucleus
    if(nucleus->getStore()->getBook()->getCascading()==0 &&
        nucleus->getStore()->getIncomingParticles().empty()) {
      DEBUG("No participants in the nucleus and no incoming particles left, stopping cascade" << std::endl);
      return false;
    }
    // Stop if the remnant is smaller than minRemnantSize
    if(nucleus->getA() <= minRemnantSize) {
      DEBUG("Remnant size (" << nucleus->getA()
          << ") smaller than or equal to minimum (" << minRemnantSize
          << "), stopping cascade" << std::endl);
      return false;
    }
    // Stop if we have to try and make a compound nucleus or if we have to
    // force a transparent
    if(nucleus->getTryCompoundNucleus()) {
      DEBUG("Trying to make a compound nucleus, stopping cascade" << std::endl);
      return false;
    }
    if(nucleus->isForcedTransparent()) {
      DEBUG("Forcing a transparent, stopping cascade" << std::endl);
      return false;
    }

    return true;
  }

  void INCL::finalizeGlobalInfo() {
    theGlobalInfo.nucleonAbsorptionCrossSection = theGlobalInfo.geometricCrossSection *
      ((G4double) theGlobalInfo.nNucleonAbsorptions) / ((G4double) theGlobalInfo.nShots);
    theGlobalInfo.pionAbsorptionCrossSection = theGlobalInfo.geometricCrossSection *
      ((G4double) theGlobalInfo.nPionAbsorptions) / ((G4double) theGlobalInfo.nShots);
    theGlobalInfo.reactionCrossSection = theGlobalInfo.geometricCrossSection *
      ((G4double) (theGlobalInfo.nShots - theGlobalInfo.nTransparents)) /
      ((G4double) theGlobalInfo.nShots);
    theGlobalInfo.errorReactionCrossSection = theGlobalInfo.geometricCrossSection *
      std::sqrt((G4double) (theGlobalInfo.nShots - theGlobalInfo.nTransparents)) /
      ((G4double) theGlobalInfo.nShots);
  }

  G4int INCL::makeProjectileRemnant() {
    G4int nUnmergedSpectators = 0;

    // Do nothing if this is not a nucleus-nucleus reaction
    if(!nucleus->getProjectileRemnant())
      return 0;

    // Get the spectators (geometrical+dynamical) from the Store
    ParticleList geomSpectators(nucleus->getProjectileRemnant()->getParticles());
    ParticleList dynSpectators(nucleus->getStore()->extractDynamicalSpectators());

    // If there are no spectators, do nothing
    if(dynSpectators.empty() && geomSpectators.empty()) {
      nucleus->deleteProjectileRemnant();
      return 0;
    } else if(geomSpectators.size()+dynSpectators.size()==1) {
      if(dynSpectators.empty()) {
        // No dynamical spectators, one geometrical spectator
        // It should already be on shell.
#if !defined(NDEBUG) && !defined(INCLXX_IN_GEANT4_MODE)
        Particle *theSpectator = geomSpectators.front();
#endif
// assert(std::abs(theSpectator->getTableMass()-theSpectator->getInvariantMass())<1.e-3);
        nucleus->moveProjectileRemnantComponentsToOutgoing();
      } else {
        // No geometrical spectators, one dynamical spectator
        // Just put it back in the outgoing list
        nucleus->getStore()->addToOutgoing(dynSpectators.front());
      }
      nucleus->deleteProjectileRemnant();
    } else {
      // Make a cluster out of the geometrical spectators
      ProjectileRemnant *theProjectileRemnant = nucleus->getProjectileRemnant();

      // Add the dynamical spectators to the bunch
      ParticleList rejected = theProjectileRemnant->addMostDynamicalSpectators(dynSpectators);
      // Put back the rejected spectators into the outgoing list
      nUnmergedSpectators = rejected.size();
      nucleus->getStore()->addToOutgoing(rejected);

      // Deal with the projectile remnant
      nucleus->finalizeProjectileRemnant(propagationModel->getCurrentTime());

    }

    return nUnmergedSpectators;
  }

  void INCL::initMaxInteractionDistance(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy) {
    if(projectileSpecies.theType != Composite) {
      maxInteractionDistance = 0.;
      return;
    }

    const G4double projectileKineticEnergyPerNucleon = kineticEnergy/projectileSpecies.theA;
    const G4double r0 = NuclearDensityFactory::createDensity(theA, theZ)->getNuclearRadius();

    maxInteractionDistance = r0 + CrossSections::interactionDistanceNN(projectileKineticEnergyPerNucleon);
  }

  void INCL::initUniverseRadius(ParticleSpecies const &p, const G4double kineticEnergy, const G4int A, const G4int Z) {
    G4double rMax = 0.0;
    if(A==0) {
      IsotopicDistribution const &anIsotopicDistribution =
        ParticleTable::getNaturalIsotopicDistribution(Z);
      IsotopeVector theIsotopes = anIsotopicDistribution.getIsotopes();
      for(IsotopeIter i=theIsotopes.begin(); i!=theIsotopes.end(); ++i) {
        NuclearDensity *theDensity = NuclearDensityFactory::createDensity(i->theA,Z);
        if(!theDensity) {
          FATAL("NULL density in initUniverseRadius. "
                << "Projectile type=" << p.theType
                << ", A=" << p.theA
                << ", Z=" << p.theZ
                << ", kinE=" << kineticEnergy
                << ", target A=" << A
                << ", Z=" << Z
                << std::endl);
        }
        rMax = std::max(theDensity->getMaximumRadius(), rMax);
      }
    } else {
      NuclearDensity *theDensity = NuclearDensityFactory::createDensity(A,Z);
      if(!theDensity) {
        FATAL("NULL density in initUniverseRadius. "
              << "Projectile type=" << p.theType
              << ", A=" << p.theA
              << ", Z=" << p.theZ
              << ", kinE=" << kineticEnergy
              << ", target A=" << A
              << ", Z=" << Z
              << std::endl);
      }
      rMax = theDensity->getMaximumRadius();
    }
    if(p.theType==Composite) {
      maxUniverseRadius = rMax;
    } else if(p.theType==Proton || p.theType==Neutron) {
      const G4double interactionDistanceNN = CrossSections::interactionDistanceNN(kineticEnergy);
      if(interactionDistanceNN>CrossSections::interactionDistanceNN1GeV()) {
        maxUniverseRadius = rMax
          - CrossSections::interactionDistanceNN1GeV()
          + interactionDistanceNN;
      } else
        maxUniverseRadius = rMax;
    } else if(p.theType==PiPlus
        || p.theType==PiZero
        || p.theType==PiMinus) {
      const G4double interactionDistancePiN = CrossSections::interactionDistancePiN(kineticEnergy);
      if(interactionDistancePiN>CrossSections::interactionDistancePiN1GeV()) {
        maxUniverseRadius = rMax
          - CrossSections::interactionDistancePiN1GeV()
          + interactionDistancePiN;
      } else
        maxUniverseRadius = rMax;
    }
  }

}
