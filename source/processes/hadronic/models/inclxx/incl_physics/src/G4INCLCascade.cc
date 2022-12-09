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

/** \file G4INCLCascade.cc
 *
 * INCL Cascade
 */
#include "G4INCLCascade.hh"
#include "G4INCLRandom.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNuclearMassTable.hh"
#include "G4INCLGlobalInfo.hh"

#include "G4INCLPauliBlocking.hh"

#include "G4INCLCrossSections.hh"

#include "G4INCLPhaseSpaceGenerator.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"

#include "G4INCLINuclearPotential.hh"

#include "G4INCLCoulombDistortion.hh"

#include "G4INCLClustering.hh"

#include "G4INCLIntersection.hh"

#include "G4INCLBinaryCollisionAvatar.hh"

#include "G4INCLCascadeAction.hh"
#include "G4INCLAvatarDumpAction.hh"

#include <cstring> 
#include <cstdlib>
#include <numeric>

namespace G4INCL {
  
  INCL::INCL(Config const * const config)
    :propagationModel(0), theA(208), theZ(82), theS(0),
    targetInitSuccess(false),
    maxImpactParameter(0.),
    maxUniverseRadius(0.),
    maxInteractionDistance(0.),
    fixedImpactParameter(0.),
    theConfig(config),
    nucleus(NULL),
    forceTransparent(false),
    minRemnantSize(4)
  {
    // Set the logger object.
#ifdef INCLXX_IN_GEANT4_MODE
    Logger::initVerbosityLevelFromEnvvar();
#else // INCLXX_IN_GEANT4_MODE
    Logger::initialize(theConfig);
#endif // INCLXX_IN_GEANT4_MODE

    // Set the random number generator algorithm. The system can support
    // multiple different generator algorithms in a completely
    // transparent way.
    Random::initialize(theConfig);

    // Select the Pauli and CDPP blocking algorithms
    Pauli::initialize(theConfig);

    // Set the cross-section set
    CrossSections::initialize(theConfig);

    // Set the phase-space generator
    PhaseSpaceGenerator::initialize(theConfig);

    // Select the Coulomb-distortion algorithm:
    CoulombDistortion::initialize(theConfig);

    // Select the clustering algorithm:
    Clustering::initialize(theConfig);

    // Initialize the INCL particle table:
    ParticleTable::initialize(theConfig);

    // Initialize the value of cutNN in BinaryCollisionAvatar
    BinaryCollisionAvatar::setCutNN(theConfig->getCutNN());

    // Initialize the value of strange cross section bias
    BinaryCollisionAvatar::setBias(theConfig->getBias());

    // Propagation model is responsible for finding avatars and
    // transporting the particles. In principle this step is "hidden"
    // behind an abstract interface and the rest of the system does not
    // care how the transportation and avatar finding is done. This
    // should allow us to "easily" experiment with different avatar
    // finding schemes and even to support things like curved
    // trajectories in the future.
    propagationModel = new StandardPropagationModel(theConfig->getLocalEnergyBBType(),theConfig->getLocalEnergyPiType(),theConfig->getHadronizationTime());
    if(theConfig->getCascadeActionType() == AvatarDumpActionType)
      cascadeAction = new AvatarDumpAction();
    else
      cascadeAction = new CascadeAction();
    cascadeAction->beforeRunAction(theConfig);

    theGlobalInfo.cascadeModel = theConfig->getVersionString();
    theGlobalInfo.deexcitationModel = theConfig->getDeExcitationString();
#ifdef INCL_ROOT_USE
    theGlobalInfo.rootSelection = theConfig->getROOTSelectionString();
#endif

#ifndef INCLXX_IN_GEANT4_MODE
    // Fill in the global information
    theGlobalInfo.At = theConfig->getTargetA();
    theGlobalInfo.Zt = theConfig->getTargetZ();
    theGlobalInfo.St = theConfig->getTargetS();
    const ParticleSpecies theSpecies = theConfig->getProjectileSpecies();
    theGlobalInfo.Ap = theSpecies.theA;
    theGlobalInfo.Zp = theSpecies.theZ;
    theGlobalInfo.Sp = theSpecies.theS;
    theGlobalInfo.Ep = theConfig->getProjectileKineticEnergy();
    theGlobalInfo.biasFactor = theConfig->getBias();
#endif

    fixedImpactParameter = theConfig->getImpactParameter();
  }

  INCL::~INCL() {
    InteractionAvatar::deleteBackupParticles();
#ifndef INCLXX_IN_GEANT4_MODE
    NuclearMassTable::deleteTable();
#endif
    PhaseSpaceGenerator::deletePhaseSpaceGenerator();
    CrossSections::deleteCrossSections();
    Pauli::deleteBlockers();
    CoulombDistortion::deleteCoulomb();
    Random::deleteGenerator();
    Clustering::deleteClusteringModel();
#ifndef INCLXX_IN_GEANT4_MODE
    Logger::deleteLoggerSlave();
#endif
    NuclearDensityFactory::clearCache();
    NuclearPotential::clearCache();
    cascadeAction->afterRunAction();
    delete cascadeAction;
    delete propagationModel;
    delete theConfig;
  }

  G4bool INCL::prepareReaction(const ParticleSpecies &projectileSpecies, const G4double kineticEnergy, const G4int A, const G4int Z, const G4int S) {
    if(A < 0 || A > 300 || Z < 1 || Z > 200) {
      INCL_ERROR("Unsupported target: A = " << A << " Z = " << Z << " S = " << S << '\n'
                 << "Target configuration rejected." << '\n');
      return false;
    }
    if(projectileSpecies.theType==Composite &&
       (projectileSpecies.theZ==projectileSpecies.theA || projectileSpecies.theZ==0)) {
      INCL_ERROR("Unsupported projectile: A = " << projectileSpecies.theA << " Z = " << projectileSpecies.theZ << " S = " << projectileSpecies.theS << '\n'
                 << "Projectile configuration rejected." << '\n');
      return false;
    }

    // Reset the forced-transparent flag
    forceTransparent = false;

    // Initialise the maximum universe radius
    initUniverseRadius(projectileSpecies, kineticEnergy, A, Z);

    // Initialise the nucleus
    theZ = Z;
    theS = S;
    if(theConfig->isNaturalTarget())
      theA = ParticleTable::drawRandomNaturalIsotope(Z);
    else
      theA = A;
    initializeTarget(theA, theZ, theS);

    // Set the maximum impact parameter
    maxImpactParameter = CoulombDistortion::maxImpactParameter(projectileSpecies, kineticEnergy, nucleus);
    INCL_DEBUG("Maximum impact parameter initialised: " << maxImpactParameter << '\n');

    // For forced CN events
    initMaxInteractionDistance(projectileSpecies, kineticEnergy);

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

  G4bool INCL::initializeTarget(const G4int A, const G4int Z, const G4int S) {
    delete nucleus;

    nucleus = new Nucleus(A, Z, S, theConfig, maxUniverseRadius);
    nucleus->getStore()->getBook().reset();
    nucleus->initializeParticles();

    propagationModel->setNucleus(nucleus);
    return true;
  }

  const EventInfo &INCL::processEvent(
      ParticleSpecies const &projectileSpecies,
      const G4double kineticEnergy,
      const G4int targetA,
      const G4int targetZ,
      const G4int targetS
      ) {
    // ReInitialize the bias vector
    Particle::INCLBiasVector.clear();
    //Particle::INCLBiasVector.Clear();
    Particle::nextBiasedCollisionID = 0;
    
    // Set the target and the projectile
    targetInitSuccess = prepareReaction(projectileSpecies, kineticEnergy, targetA, targetZ, targetS);

    if(!targetInitSuccess) {
      INCL_WARN("Target initialisation failed for A=" << targetA << ", Z=" << targetZ << ", S=" << targetS << '\n');
      theEventInfo.transparent=true;
      return theEventInfo;
    }

    cascadeAction->beforeCascadeAction(propagationModel);

    const G4bool canRunCascade = preCascade(projectileSpecies, kineticEnergy);
    if(canRunCascade) {
      cascade();
      postCascade();
      cascadeAction->afterCascadeAction(nucleus);
    }
    updateGlobalInfo();
    return theEventInfo;
  }

  G4bool INCL::preCascade(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy) {
    // Reset theEventInfo
    theEventInfo.reset();

    EventInfo::eventNumber++;

    // Fill in the event information
    theEventInfo.projectileType = projectileSpecies.theType;
    theEventInfo.Ap = (G4INCL::Short_t)projectileSpecies.theA;
    theEventInfo.Zp = (G4INCL::Short_t)projectileSpecies.theZ;
    theEventInfo.Sp = (G4INCL::Short_t)projectileSpecies.theS;
    theEventInfo.Ep = kineticEnergy;
    theEventInfo.At = (G4INCL::Short_t)nucleus->getA();
    theEventInfo.Zt = (G4INCL::Short_t)nucleus->getZ();
    theEventInfo.St = (G4INCL::Short_t)nucleus->getS();

    // Do nothing below the Coulomb barrier
    if(maxImpactParameter<=0.) {
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
    INCL_DEBUG("Selected impact parameter: " << impactParameter << '\n');

    // Fill in the event information
    theEventInfo.impactParameter = impactParameter;

    const G4double effectiveImpactParameter = propagationModel->shoot(projectileSpecies, kineticEnergy, impactParameter, phi);
    if(effectiveImpactParameter < 0.) {
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
    FinalState *finalState = new FinalState;

    unsigned long loopCounter = 0;
    const unsigned long maxLoopCounter = 10000000;
    do {
      // Run book keeping actions that should take place before propagation:
      cascadeAction->beforePropagationAction(propagationModel);

      // Get the avatar with the smallest time and propagate particles
      // to that point in time.
      IAvatar *avatar = propagationModel->propagate(finalState);

      finalState->reset();

      // Run book keeping actions that should take place after propagation:
      cascadeAction->afterPropagationAction(propagationModel, avatar);

      if(avatar == 0) break; // No more avatars in the avatar list.

      // Run book keeping actions that should take place before avatar:
      cascadeAction->beforeAvatarAction(avatar, nucleus);

      // Channel is responsible for calculating the outcome of the
      // selected avatar. There are different kinds of channels. The
      // class IChannel is, again, an abstract interface that defines
      // the externally observable behavior of all interaction
      // channels.
      // The handling of the channel is transparent to the API.
      // Final state tells what changed...
      avatar->fillFinalState(finalState);
      // Run book keeping actions that should take place after avatar:
      cascadeAction->afterAvatarAction(avatar, nucleus, finalState);

      // So now we must give this information to the nucleus
      nucleus->applyFinalState(finalState);
      // and now we are ready to process the next avatar!

      delete avatar;

      ++loopCounter;
    } while(continueCascade() && loopCounter<maxLoopCounter); /* Loop checking, 10.07.2015, D.Mancusi */
    
    delete finalState;
  }

  void INCL::postCascade() {
    // Fill in the event information
    theEventInfo.stoppingTime = propagationModel->getCurrentTime();

    // The event bias
    theEventInfo.eventBias = (Double_t) Particle::getTotalBias();
    
    // Forced CN?
    if(nucleus->getTryCompoundNucleus()) {
      INCL_DEBUG("Trying compound nucleus" << '\n');
      makeCompoundNucleus();
      theEventInfo.transparent = forceTransparent;
      // Global checks of conservation laws
#ifndef INCLXX_IN_GEANT4_MODE
      if(!theEventInfo.transparent) globalConservationChecks(true);
#endif
      return;
    }

    theEventInfo.transparent = forceTransparent || nucleus->isEventTransparent();

    if(theEventInfo.transparent) {
      ProjectileRemnant * const projectileRemnant = nucleus->getProjectileRemnant();
      if(projectileRemnant) {
        // Clear the incoming list (particles will be deleted by the ProjectileRemnant)
        nucleus->getStore()->clearIncoming();
      } else {
        // Delete particles in the incoming list
        nucleus->getStore()->deleteIncoming();
      }
    } else {
      
      // Check if the nucleus contains strange particles
      theEventInfo.sigmasInside = nucleus->containsSigma();
      theEventInfo.antikaonsInside = nucleus->containsAntiKaon();
      theEventInfo.lambdasInside = nucleus->containsLambda();
      theEventInfo.kaonsInside = nucleus->containsKaon();
      
      // Capture antiKaons and Sigmas and produce Lambda instead
      theEventInfo.absorbedStrangeParticle = nucleus->decayInsideStrangeParticles();
      
      // Emit strange particles still inside the nucleus
      nucleus->emitInsideStrangeParticles();
      theEventInfo.emitKaon = nucleus->emitInsideKaon();

#ifdef INCLXX_IN_GEANT4_MODE
      theEventInfo.emitLambda = nucleus->emitInsideLambda();
#endif // INCLXX_IN_GEANT4_MODE
      
      // Check if the nucleus contains deltas
      theEventInfo.deltasInside = nucleus->containsDeltas();

      // Take care of any remaining deltas
      theEventInfo.forcedDeltasOutside = nucleus->decayOutgoingDeltas();
      theEventInfo.forcedDeltasInside = nucleus->decayInsideDeltas();

      // Take care of any remaining etas, omegas, neutral Sigmas and/or neutral kaons
      G4double timeThreshold=theConfig->getDecayTimeThreshold();
      theEventInfo.forcedPionResonancesOutside = nucleus->decayOutgoingPionResonances(timeThreshold);
      nucleus->decayOutgoingSigmaZero(timeThreshold);
      nucleus->decayOutgoingNeutralKaon();
        
      // Apply Coulomb distortion, if appropriate
      // Note that this will apply Coulomb distortion also on pions emitted by
      // unphysical remnants (see decayInsideDeltas). This is at variance with
      // what INCL4.6 does, but these events are (should be!) so rare that
      // whatever we do doesn't (shouldn't!) make any noticeable difference.
      CoulombDistortion::distortOut(nucleus->getStore()->getOutgoingParticles(), nucleus);

      // If the normal cascade predicted complete fusion, use the tabulated
      // masses to compute the excitation energy, the recoil, etc.
      if(nucleus->getStore()->getOutgoingParticles().size()==0
         && (!nucleus->getProjectileRemnant()
             || nucleus->getProjectileRemnant()->getParticles().size()==0)) {

        INCL_DEBUG("Cascade resulted in complete fusion, using realistic fusion kinematics" << '\n');

        nucleus->useFusionKinematics();

        if(nucleus->getExcitationEnergy()<0.) {
          // Complete fusion is energetically impossible, return a transparent
          INCL_WARN("Complete-fusion kinematics yields negative excitation energy, returning a transparent!" << '\n');
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
        if(nucleus->getA()==1 && minRemnantSize>1) {
          INCL_ERROR("Computing one-nucleon recoil kinematics. We should never be here nowadays, cascade should stop earlier than this." << '\n');
        }
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
      theEventInfo.clusterDecay = nucleus->decayOutgoingClusters() || nucleus->decayMe();

#ifndef INCLXX_IN_GEANT4_MODE
      // Global checks of conservation laws
      globalConservationChecks(true);
#endif

      // Fill the EventInfo structure
      nucleus->fillEventInfo(&theEventInfo);

    }
  }

  void INCL::makeCompoundNucleus() {
    // If this is not a nucleus-nucleus collision, don't attempt to make a
    // compound nucleus.
    //
    // Yes, even nucleon-nucleus collisions can lead to particles entering
    // below the Fermi level. Take e.g. 1-MeV p + He4.
    if(!nucleus->isNucleusNucleusCollision()) {
      forceTransparent = true;
      return;
    }

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
    const G4double theTargetMass = ParticleTable::getTableMass(theEventInfo.At, theEventInfo.Zt, theEventInfo.St);
    G4int theCNA=theEventInfo.At, theCNZ=theEventInfo.Zt, theCNS=theEventInfo.St;
    Cluster * const theProjectileRemnant = nucleus->getProjectileRemnant();
    G4double theCNEnergy = theTargetMass + theProjectileRemnant->getEnergy();

    // Loop over the potential participants
    ParticleList const &initialProjectileComponents = theProjectileRemnant->getParticles();
    std::vector<Particle *> shuffledComponents(initialProjectileComponents.begin(), initialProjectileComponents.end());
    // Shuffle the list of potential participants
    std::shuffle(shuffledComponents.begin(), shuffledComponents.end(), Random::getAdapter());

    G4bool success = true;
    G4bool atLeastOneNucleonEntering = false;
    for(std::vector<Particle*>::const_iterator p=shuffledComponents.begin(), e=shuffledComponents.end(); p!=e; ++p) {
      // Skip particles that miss the interaction distance
      Intersection intersectionInteractionDistance(IntersectionFactory::getEarlierTrajectoryIntersection(
            (*p)->getPosition(),
            (*p)->getPropagationVelocity(),
            maxInteractionDistance));
      if(!intersectionInteractionDistance.exists)
        continue;

      // Build an entry avatar for this nucleon
      atLeastOneNucleonEntering = true;
      ParticleEntryAvatar *theAvatar = new ParticleEntryAvatar(0.0, nucleus, *p);
      nucleus->getStore()->addParticleEntryAvatar(theAvatar);
      FinalState *fs = theAvatar->getFinalState();
      nucleus->applyFinalState(fs);
      FinalStateValidity validity = fs->getValidity();
      delete fs;
      switch(validity) {
        case ValidFS:
        case ParticleBelowFermiFS:
        case ParticleBelowZeroFS:
          // Add the particle to the CN
          theCNA++;
          theCNZ += (*p)->getZ();
          theCNS += (*p)->getS();
          break;
        case PauliBlockedFS:
        case NoEnergyConservationFS:
        default:
          success = false;
          break;
      }
    }

    if(!success || !atLeastOneNucleonEntering) {
      INCL_DEBUG("No nucleon entering in forced CN, forcing a transparent" << '\n');
      forceTransparent = true;
      return;
    }

// assert(theCNA==nucleus->getA());
// assert(theCNA<=theEventInfo.At+theEventInfo.Ap);
// assert(theCNZ<=theEventInfo.Zt+theEventInfo.Zp);
// assert(theCNS>=theEventInfo.St+theEventInfo.Sp);

    // Update the kinematics of the CN
    theCNEnergy -= theProjectileRemnant->getEnergy();
    theCNMomentum -= theProjectileRemnant->getMomentum();

    // Deal with the projectile remnant
    nucleus->finalizeProjectileRemnant(propagationModel->getCurrentTime());

    // Subtract the angular momentum of the projectile remnant
// assert(nucleus->getStore()->getOutgoingParticles().empty());
    theCNSpin -= theProjectileRemnant->getAngularMomentum();

    // Compute the excitation energy of the CN
    const G4double theCNMass = ParticleTable::getTableMass(theCNA,theCNZ,theCNS);
    const G4double theCNInvariantMassSquared = theCNEnergy*theCNEnergy-theCNMomentum.mag2();
    if(theCNInvariantMassSquared<0.) {
      // Negative invariant mass squared, return a transparent
      forceTransparent = true;
      return;
    }
    const G4double theCNExcitationEnergy = std::sqrt(theCNInvariantMassSquared) - theCNMass;
    if(theCNExcitationEnergy<0.) {
      // Negative excitation energy, return a transparent
      INCL_DEBUG("CN excitation energy is negative, forcing a transparent" << '\n'
            << "  theCNA = " << theCNA << '\n'
            << "  theCNZ = " << theCNZ << '\n'
            << "  theCNS = " << theCNS << '\n'
            << "  theCNEnergy = " << theCNEnergy << '\n'
            << "  theCNMomentum = (" << theCNMomentum.getX() << ", "<< theCNMomentum.getY() << ", "  << theCNMomentum.getZ() << ")" << '\n'
            << "  theCNExcitationEnergy = " << theCNExcitationEnergy << '\n'
            << "  theCNSpin = (" << theCNSpin.getX() << ", "<< theCNSpin.getY() << ", "  << theCNSpin.getZ() << ")" << '\n'
            );
      forceTransparent = true;
      return;
    } else {
      // Positive excitation energy, can make a CN
      INCL_DEBUG("CN excitation energy is positive, forcing a CN" << '\n'
            << "  theCNA = " << theCNA << '\n'
            << "  theCNZ = " << theCNZ << '\n'
            << "  theCNS = " << theCNS << '\n'
            << "  theCNEnergy = " << theCNEnergy << '\n'
            << "  theCNMomentum = (" << theCNMomentum.getX() << ", "<< theCNMomentum.getY() << ", "  << theCNMomentum.getZ() << ")" << '\n'
            << "  theCNExcitationEnergy = " << theCNExcitationEnergy << '\n'
            << "  theCNSpin = (" << theCNSpin.getX() << ", "<< theCNSpin.getY() << ", "  << theCNSpin.getZ() << ")" << '\n'
            );
      nucleus->setA(theCNA);
      nucleus->setZ(theCNZ);
      nucleus->setS(theCNS);
      nucleus->setMomentum(theCNMomentum);
      nucleus->setEnergy(theCNEnergy);
      nucleus->setExcitationEnergy(theCNExcitationEnergy);
      nucleus->setMass(theCNMass+theCNExcitationEnergy);
      nucleus->setSpin(theCNSpin); // neglects any orbital angular momentum of the CN

      // Take care of any remaining deltas
      theEventInfo.forcedDeltasOutside = nucleus->decayOutgoingDeltas();

      // Take care of any remaining etas and/or omegas
      G4double timeThreshold=theConfig->getDecayTimeThreshold();
      theEventInfo.forcedPionResonancesOutside = nucleus->decayOutgoingPionResonances(timeThreshold);
      
      // Take care of any remaining Kaons
      theEventInfo.emitKaon = nucleus->emitInsideKaon();
        
      // Cluster decay
      theEventInfo.clusterDecay = nucleus->decayOutgoingClusters() || nucleus->decayMe();

      // Fill the EventInfo structure
      nucleus->fillEventInfo(&theEventInfo);
    }
  }

  void INCL::rescaleOutgoingForRecoil() {
    RecoilCMFunctor theRecoilFunctor(nucleus, theEventInfo);

    // Apply the root-finding algorithm
    const RootFinder::Solution theSolution = RootFinder::solve(&theRecoilFunctor, 1.0);
    if(theSolution.success) {
      theRecoilFunctor(theSolution.x); // Apply the solution
    } else {
      INCL_WARN("Couldn't accommodate remnant recoil while satisfying energy conservation, root-finding algorithm failed." << '\n');
    }

  }

#ifndef INCLXX_IN_GEANT4_MODE
  void INCL::globalConservationChecks(G4bool afterRecoil) {
    Nucleus::ConservationBalance theBalance = nucleus->getConservationBalance(theEventInfo,afterRecoil);

    // Global conservation checks
    const G4double pLongBalance = theBalance.momentum.getZ();
    const G4double pTransBalance = theBalance.momentum.perp();
    if(theBalance.Z != 0) {
      INCL_ERROR("Violation of charge conservation! ZBalance = " << theBalance.Z << " eventNumber=" << theEventInfo.eventNumber << '\n');
    }
    if(theBalance.A != 0) {
      INCL_ERROR("Violation of baryon-number conservation! ABalance = " << theBalance.A << " Emit Lambda=" << theEventInfo.emitLambda << " eventNumber=" << theEventInfo.eventNumber << '\n');
    }
    if(theBalance.S != 0) {
      INCL_ERROR("Violation of strange-number conservation! SBalance = " << theBalance.S << " eventNumber=" << theEventInfo.eventNumber << '\n');
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
      INCL_WARN("Violation of energy conservation > " << EThreshold << " MeV. EBalance = " << theBalance.energy << " Emit Lambda=" << theEventInfo.emitLambda << " afterRecoil = " << afterRecoil << " eventNumber=" << theEventInfo.eventNumber << '\n');
    }
    if(std::abs(pLongBalance)>pLongThreshold) {
      INCL_WARN("Violation of longitudinal momentum conservation > " << pLongThreshold << " MeV/c. pLongBalance = " << pLongBalance << " afterRecoil = " << afterRecoil << " eventNumber=" << theEventInfo.eventNumber << '\n');
    }
    if(std::abs(pTransBalance)>pTransThreshold) {
      INCL_WARN("Violation of transverse momentum conservation > " << pTransThreshold << " MeV/c. pTransBalance = " << pTransBalance << " afterRecoil = " << afterRecoil << " eventNumber=" << theEventInfo.eventNumber << '\n');
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
      INCL_DEBUG("Cascade time (" << propagationModel->getCurrentTime()
          << ") exceeded stopping time (" << propagationModel->getStoppingTime()
          << "), stopping cascade" << '\n');
      return false;
    }
    // Stop if there are no participants and no pions inside the nucleus
    if(nucleus->getStore()->getBook().getCascading()==0 &&
        nucleus->getStore()->getIncomingParticles().empty()) {
      INCL_DEBUG("No participants in the nucleus and no incoming particles left, stopping cascade" << '\n');
      return false;
    }
    // Stop if the remnant is smaller than minRemnantSize
    if(nucleus->getA() <= minRemnantSize) {
      INCL_DEBUG("Remnant size (" << nucleus->getA()
          << ") smaller than or equal to minimum (" << minRemnantSize
          << "), stopping cascade" << '\n');
      return false;
    }
    // Stop if we have to try and make a compound nucleus or if we have to
    // force a transparent
    if(nucleus->getTryCompoundNucleus()) {
      INCL_DEBUG("Trying to make a compound nucleus, stopping cascade" << '\n');
      return false;
    }

    return true;
  }

  void INCL::finalizeGlobalInfo(Random::SeedVector const &initialSeeds) {
    const G4double normalisationFactor = theGlobalInfo.geometricCrossSection /
      ((G4double) theGlobalInfo.nShots);
    theGlobalInfo.nucleonAbsorptionCrossSection = normalisationFactor *
      ((G4double) theGlobalInfo.nNucleonAbsorptions);
    theGlobalInfo.pionAbsorptionCrossSection = normalisationFactor *
      ((G4double) theGlobalInfo.nPionAbsorptions);
    theGlobalInfo.reactionCrossSection = normalisationFactor *
      ((G4double) (theGlobalInfo.nShots - theGlobalInfo.nTransparents));
    theGlobalInfo.errorReactionCrossSection = normalisationFactor *
      std::sqrt((G4double) (theGlobalInfo.nShots - theGlobalInfo.nTransparents));
    theGlobalInfo.forcedCNCrossSection = normalisationFactor *
      ((G4double) theGlobalInfo.nForcedCompoundNucleus);
    theGlobalInfo.errorForcedCNCrossSection = normalisationFactor *
      std::sqrt((G4double) (theGlobalInfo.nForcedCompoundNucleus));
    theGlobalInfo.completeFusionCrossSection = normalisationFactor *
      ((G4double) theGlobalInfo.nCompleteFusion);
    theGlobalInfo.errorCompleteFusionCrossSection = normalisationFactor *
      std::sqrt((G4double) (theGlobalInfo.nCompleteFusion));
    theGlobalInfo.energyViolationInteractionCrossSection = normalisationFactor *
      ((G4double) theGlobalInfo.nEnergyViolationInteraction);

    theGlobalInfo.initialRandomSeeds.assign(initialSeeds.begin(), initialSeeds.end());

    Random::SeedVector theSeeds = Random::getSeeds();
    theGlobalInfo.finalRandomSeeds.assign(theSeeds.begin(), theSeeds.end());
  }

  G4int INCL::makeProjectileRemnant() {
    // Do nothing if this is not a nucleus-nucleus reaction
    if(!nucleus->getProjectileRemnant())
      return 0;

    // Get the spectators (geometrical+dynamical) from the Store
    ParticleList geomSpectators(nucleus->getProjectileRemnant()->getParticles());
    ParticleList dynSpectators(nucleus->getStore()->extractDynamicalSpectators());

    G4int nUnmergedSpectators = 0;

    // If there are no spectators, do nothing
    if(dynSpectators.empty() && geomSpectators.empty()) {
      return 0;
    } else if(dynSpectators.size()==1 && geomSpectators.empty()) {
      // No geometrical spectators, one dynamical spectator
      // Just put it back in the outgoing list
      nucleus->getStore()->addToOutgoing(dynSpectators.front());
    } else {
      // Make a cluster out of the geometrical spectators
      ProjectileRemnant *theProjectileRemnant = nucleus->getProjectileRemnant();

      // Add the dynamical spectators to the bunch
      ParticleList rejected = theProjectileRemnant->addAllDynamicalSpectators(dynSpectators);
      // Put back the rejected spectators into the outgoing list
      nUnmergedSpectators = (G4int)rejected.size();
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

    const G4double r0 = std::max(ParticleTable::getNuclearRadius(Proton, theA, theZ),
                               ParticleTable::getNuclearRadius(Neutron, theA, theZ));

    const G4double theNNDistance = CrossSections::interactionDistanceNN(projectileSpecies, kineticEnergy);
    maxInteractionDistance = r0 + theNNDistance;
    INCL_DEBUG("Initialised interaction distance: r0 = " << r0 << '\n'
          << "    theNNDistance = " << theNNDistance << '\n'
          << "    maxInteractionDistance = " << maxInteractionDistance << '\n');
  }

  void INCL::initUniverseRadius(ParticleSpecies const &p, const G4double kineticEnergy, const G4int A, const G4int Z) {
    G4double rMax = 0.0;
    if(A==0) {
      IsotopicDistribution const &anIsotopicDistribution =
        ParticleTable::getNaturalIsotopicDistribution(Z);
      IsotopeVector theIsotopes = anIsotopicDistribution.getIsotopes();
      for(IsotopeIter i=theIsotopes.begin(), e=theIsotopes.end(); i!=e; ++i) {
        const G4double pMaximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, i->theA, Z);
        const G4double nMaximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, i->theA, Z);
        const G4double maximumRadius = std::max(pMaximumRadius, nMaximumRadius);
        rMax = std::max(maximumRadius, rMax);
      }
    } else {
      const G4double pMaximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, A, Z);
      const G4double nMaximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, A, Z);
      const G4double maximumRadius = std::max(pMaximumRadius, nMaximumRadius);
      rMax = std::max(maximumRadius, rMax);
    }
    if(p.theType==Composite || p.theType==Proton || p.theType==Neutron) {
      const G4double interactionDistanceNN = CrossSections::interactionDistanceNN(p, kineticEnergy);
      maxUniverseRadius = rMax + interactionDistanceNN;
    } else if(p.theType==PiPlus
        || p.theType==PiZero
        || p.theType==PiMinus) {
      const G4double interactionDistancePiN = CrossSections::interactionDistancePiN(kineticEnergy);
      maxUniverseRadius = rMax + interactionDistancePiN;
    } else if(p.theType==KPlus
        || p.theType==KZero) {
      const G4double interactionDistanceKN = CrossSections::interactionDistanceKN(kineticEnergy);
      maxUniverseRadius = rMax + interactionDistanceKN;
    } else if(p.theType==KZeroBar
        || p.theType==KMinus) {
      const G4double interactionDistanceKbarN = CrossSections::interactionDistanceKbarN(kineticEnergy);
      maxUniverseRadius = rMax + interactionDistanceKbarN;
    } else if(p.theType==Lambda
        ||p.theType==SigmaPlus
        || p.theType==SigmaZero
        || p.theType==SigmaMinus) {
      const G4double interactionDistanceYN = CrossSections::interactionDistanceYN(kineticEnergy);
      maxUniverseRadius = rMax + interactionDistanceYN;
    }
    INCL_DEBUG("Initialised universe radius: " << maxUniverseRadius << '\n');
  }

  void INCL::updateGlobalInfo() {
    // Increment the global counter for the number of shots
    theGlobalInfo.nShots++;

    if(theEventInfo.transparent) {
      // Increment the global counter for the number of transparents
      theGlobalInfo.nTransparents++;
      // Increment the global counter for the number of forced transparents
      if(forceTransparent)
        theGlobalInfo.nForcedTransparents++;
      return;
    }

    // Check if we have an absorption:
    if(theEventInfo.nucleonAbsorption) theGlobalInfo.nNucleonAbsorptions++;
    if(theEventInfo.pionAbsorption) theGlobalInfo.nPionAbsorptions++;

    // Count complete-fusion events
    if(theEventInfo.nCascadeParticles==0) theGlobalInfo.nCompleteFusion++;

    if(nucleus->getTryCompoundNucleus())
      theGlobalInfo.nForcedCompoundNucleus++;

    // Counters for the number of violations of energy conservation in
    // collisions
    theGlobalInfo.nEnergyViolationInteraction += theEventInfo.nEnergyViolationInteraction;
  }

}
