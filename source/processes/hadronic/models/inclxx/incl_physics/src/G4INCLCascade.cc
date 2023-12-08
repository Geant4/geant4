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
#include "G4INCLNucleus.hh"

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

#include "G4INCLPbarAtrestEntryChannel.hh"

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

//D
    //reset
    G4bool ProtonIsTheVictim = false; 
    G4bool NeutronIsTheVictim = false;
    theEventInfo.annihilationP = false;
    theEventInfo.annihilationN = false;

    //G4double AnnihilationBarrier = kineticEnergy;
    if(projectileSpecies.theType == antiProton && kineticEnergy <= theConfig->getAtrestThreshold()){
      G4double SpOverSn = 1.331;//from experiments with deuteron (E.Klempt)
      //INCL_WARN("theA number set to A-1 from " << A <<'\n');

      G4double neutronprob;
      if(theConfig->isNaturalTarget()){ // A = 0 in this case
        theA = ParticleTable::drawRandomNaturalIsotope(Z) - 1; //43 and 61 are ok (Technetium and Promethium)
        neutronprob = (theA + 1 - Z)/(theA + 1 - Z + SpOverSn*Z);  
      }
      else{
        theA = A - 1;
        neutronprob = (A - Z)/(A - Z + SpOverSn*Z);  //from experiments with deuteron (E.Klempt)
      }

      theS = S;
      
      G4double rndm = Random::shoot();
      if(rndm >= neutronprob){     //proton is annihilated
        theEventInfo.annihilationP = true;
        theZ = Z - 1;
        ProtonIsTheVictim = true;
        //INCL_WARN("theZ number set to Z-1 from " << Z << '\n');
      }  
      else{        //neutron is annihilated
        theEventInfo.annihilationN = true;
        theZ = Z;
        NeutronIsTheVictim = true;
      }  
    }
    else{ // not annihilation of pbar
      theZ = Z;
      theS = S;
      if(theConfig->isNaturalTarget())
        theA = ParticleTable::drawRandomNaturalIsotope(Z); //change order
      else
        theA = A;
    }

    AnnihilationType theAType = Def;
    if(ProtonIsTheVictim == true && NeutronIsTheVictim == false)
    theAType = PType;
    if(NeutronIsTheVictim == true && ProtonIsTheVictim == false)
    theAType = NType;

//D

    initializeTarget(theA, theZ, theS, theAType);
    
    // Set the maximum impact parameter
    maxImpactParameter = CoulombDistortion::maxImpactParameter(projectileSpecies, kineticEnergy, nucleus);
    INCL_DEBUG("Maximum impact parameter initialised: " << maxImpactParameter << '\n');

    // For forced CN events
    initMaxInteractionDistance(projectileSpecies, kineticEnergy);
// Set the geometric cross sectiony section
    if(projectileSpecies.theType == antiProton && kineticEnergy <= theConfig->getAtrestThreshold()){
      G4int currentA = A;
      if(theConfig->isNaturalTarget()){
        currentA = ParticleTable::drawRandomNaturalIsotope(Z);
      }
      G4double kineticEnergy2=kineticEnergy;
      if (kineticEnergy2 <= 0.) kineticEnergy2=0.001;
      theGlobalInfo.geometricCrossSection = 9.7* //normalization factor from Corradini
        Math::pi*std::pow((1.840 + 1.120*std::pow(currentA,(1./3.))),2)*
        (1. + (Z*G4INCL::PhysicalConstants::eSquared*(currentA+1))/(currentA*kineticEnergy2*(1.840 + 1.120*std::pow(currentA,(1./3.))))); 
         //xsection formula was borrowed from Corradini et al. https://doi.org/10.1016/j.physletb.2011.09.069    
    }
    else{
      theGlobalInfo.geometricCrossSection =
        Math::tenPi*std::pow(maxImpactParameter,2);
    }

    // Set the minimum remnant size
    if(projectileSpecies.theA > 0)
      minRemnantSize = std::min(theA, 4);
    else
      minRemnantSize = std::min(theA-1, 4);
    return true;
  }

  G4bool INCL::initializeTarget(const G4int A, const G4int Z, const G4int S, AnnihilationType theAType) { 
    delete nucleus;

    if (theAType==PType || theAType==NType) {
      G4double newmaxUniverseRadius=0.;
      if (theAType==PType) newmaxUniverseRadius=initUniverseRadiusForAntiprotonAtRest(A+1, Z+1);
      else newmaxUniverseRadius=initUniverseRadiusForAntiprotonAtRest(A+1, Z);
      nucleus = new Nucleus(A, Z, S, theConfig, newmaxUniverseRadius, theAType);
    }
    else{
      nucleus = new Nucleus(A, Z, S, theConfig, maxUniverseRadius, theAType);
    }
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

    ParticleList starlistH2;

    if (projectileSpecies.theType==antiProton && (targetA==1 || targetA==2) && targetZ==1 && targetS==0) {

      if (targetA==1) {
        preCascade_pbarH1(projectileSpecies, kineticEnergy);
      } else {
        preCascade_pbarH2(projectileSpecies, kineticEnergy);
        theEventInfo.annihilationP = false;
        theEventInfo.annihilationN = false;

        G4double SpOverSn = 1.331;  //from experiments with deuteron (E.Klempt)

        ThreeVector dummy(0.,0.,0.);
        G4double rndm = Random::shoot()*(SpOverSn+1);
        if (rndm <= SpOverSn) {  //proton is annihilated
          theEventInfo.annihilationP = true;
          Particle *p2 = new Particle(Neutron, dummy, dummy);
          starlistH2.push_back(p2);
          //delete p2;
        } else {                 //neutron is annihilated
          theEventInfo.annihilationN = true;
          Particle *p2 = new Particle(Proton, dummy, dummy);
          starlistH2.push_back(p2);
          //delete p2;
        }
      }

      // File names
#ifdef INCLXX_IN_GEANT4_MODE
      if (!G4FindDataDir("G4INCLDATA") ) {
        G4ExceptionDescription ed;
        ed << " Data missing: set environment variable G4INCLDATA\n"
           << " to point to the directory containing data files needed\n"
           << " by the INCL++ model" << G4endl;
        G4Exception("G4INCLDataFile::readData()","rawppbarFS.dat, ...", FatalException, ed);
      }
      G4String dataPath0(G4FindDataDir("G4INCLDATA"));
      G4String dataPathppbar(dataPath0 + "/rawppbarFS.dat");
      G4String dataPathnpbar(dataPath0 + "/rawnpbarFS.dat");
      G4String dataPathppbark(dataPath0 + "/rawppbarFSkaonic.dat");
      G4String dataPathnpbark(dataPath0 + "/rawnpbarFSkaonic.dat");
#else
      G4String path;
      if (theConfig) path = theConfig->getINCLXXDataFilePath();
      G4String dataPathppbar(path + "/rawppbarFS.dat");
      INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar final states" << dataPathppbar << '\n');
      G4String dataPathnpbar(path + "/rawnpbarFS.dat");
      INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar final states" << dataPathnpbar << '\n');
      G4String dataPathppbark(path + "/rawppbarFSkaonic.dat");
      INCL_DEBUG("Reading https://doi.org/10.1016/j.physrep.2005.03.002 ppbar kaonic final states" << dataPathppbark << '\n');
      G4String dataPathnpbark(path + "/rawnpbarFSkaonic.dat");
      INCL_DEBUG("Reading https://doi.org/10.1007/BF02818764 and https://link.springer.com/article/10.1007/BF02754930 npbar kaonic final states" << dataPathnpbark << '\n');
      #endif

      //read probabilities and particle types from file
      std::vector<G4double> probabilities;  //will store each FS yield
      std::vector<std::vector<G4String>> particle_types;  //will store particle names
      G4double sum = 0.0;  //will contain a sum of probabilities of all FS in the file
      G4double kaonicFSprob=0.05;  //probability to kave kaonic FS

      ParticleList starlist;
      ThreeVector mommy;  //momentum to be assigned later

      G4double rdm = Random::shoot();
      ThreeVector annihilationPosition(0.,0.,0.);
      if (rdm < (1.-kaonicFSprob)) {  // pionic FS was chosen
        INCL_DEBUG("pionic pp final state chosen" << '\n');
        sum = read_file(dataPathppbar, probabilities, particle_types);
        rdm = (rdm/(1.-kaonicFSprob))*sum;  //99.88 normalize by the sum of probabilities in the file
        //now get the line number in the file where the FS particles are stored:
        G4int n = findStringNumber(rdm, probabilities)-1;
        if ( n < 0 ) return theEventInfo;
        for (G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++) {
          if (particle_types[n][j] == "pi0") {
            Particle *p = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "pi-") {
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "pi+") {
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "omega") {
            Particle *p = new Particle(Omega, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "eta") {
            Particle *p = new Particle(Eta, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "rho-") {
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(pp);
          } else if (particle_types[n][j] == "rho+") {
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(pp);
          } else if (particle_types[n][j] == "rho0") {
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(pp);
          } else {
            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
            for (G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++) {
#ifdef INCLXX_IN_GEANT4_MODE
              G4cout << "gotcha! " << particle_types[n][jj] << G4endl;
#else
              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
#endif              
            }
#ifdef INCLXX_IN_GEANT4_MODE
            G4cout << "Some non-existing FS particle detected when reading pbar FS files" << G4endl;
#else
            std::cout << "Some non-existing FS particle detected when reading pbar FS files" << std::endl;
#endif              
          }
        }
      } else {
        INCL_DEBUG("kaonic pp final state chosen" << '\n');
        sum = read_file(dataPathppbark, probabilities, particle_types);
        rdm = ((1.-rdm)/kaonicFSprob)*sum;  //2670 normalize by the sum of probabilities in the file
        //now get the line number in the file where the FS particles are stored:
        G4int n = findStringNumber(rdm, probabilities)-1;
        if ( n < 0 ) return theEventInfo;
        for (G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++) {
          if (particle_types[n][j] == "pi0") {
            Particle *p = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "pi-") {
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "pi+") {
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "omega") {
            Particle *p = new Particle(Omega, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "eta") {
            Particle *p = new Particle(Eta, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "K-") {
            Particle *p = new Particle(KMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "K+") {
            Particle *p = new Particle(KPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "K0") {
            Particle *p = new Particle(KZero, mommy, annihilationPosition);
            starlist.push_back(p);
          } else if (particle_types[n][j] == "K0b") {
            Particle *p = new Particle(KZeroBar, mommy, annihilationPosition);
            starlist.push_back(p);
          } else {
            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
            for (G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++) {
#ifdef INCLXX_IN_GEANT4_MODE
              G4cout << "gotcha! " << particle_types[n][jj] << G4endl;
#else
              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
#endif              
            }
#ifdef INCLXX_IN_GEANT4_MODE
            G4cout << "Some non-existing FS particle detected when reading pbar FS files" << G4endl;
#else
            std::cout << "Some non-existing FS particle detected when reading pbar FS files" << std::endl;
#endif              
          }
        }
      }

      //compute energies of mesons with a phase-space model
      G4double energyOfMesonStar=ParticleTable::getRealMass(Proton)+ParticleTable::getRealMass(antiProton);
      if (starlist.size() < 2) {
        INCL_ERROR("should never happen, at least 2 final state particles!" << '\n');
      } else if (starlist.size() == 2) {
        ParticleIter first = starlist.begin();
        ParticleIter last = std::next(first, 1);
        G4double m1 = (*first)->getMass();
        G4double m2 = (*last)->getMass();
        G4double s = energyOfMesonStar*energyOfMesonStar;
        G4double mom1 = std::sqrt(s/4. - (std::pow(m1,2) + std::pow(m2,2))/2. - std::pow(m1,2)*std::pow(m2,2)/s + (std::pow(m1,4) + 2.*std::pow(m1*m2,2) + std::pow(m2,4))/(4.*s));
        ThreeVector momentello = Random::normVector(mom1);
        (*first)->setMomentum(momentello);
        (*first)->adjustEnergyFromMomentum();
        (*last)->setMomentum(-momentello);
        (*last)->adjustEnergyFromMomentum();
      } else {
        PhaseSpaceGenerator::generate(energyOfMesonStar, starlist);
      }

      if (targetA==1) postCascade_pbarH1(starlist);
      else            postCascade_pbarH2(starlist,starlistH2);

      theGlobalInfo.nShots++;
      return theEventInfo;
    }  // pbar on H1

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
      postCascade(projectileSpecies, kineticEnergy);
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
    theEventInfo.Ap = (Short_t)projectileSpecies.theA;
    theEventInfo.Zp = (Short_t)projectileSpecies.theZ;
    theEventInfo.Sp = (Short_t)projectileSpecies.theS;
    theEventInfo.Ep = kineticEnergy;
    theEventInfo.St = (Short_t)nucleus->getS();

    if(nucleus->getAnnihilationType()==PType){
      theEventInfo.annihilationP = true;
      theEventInfo.At = (Short_t)nucleus->getA()+1;
      theEventInfo.Zt = (Short_t)nucleus->getZ()+1;
    }
    else if(nucleus->getAnnihilationType()==NType){
      theEventInfo.annihilationN = true;
      theEventInfo.At = (Short_t)nucleus->getA()+1;
      theEventInfo.Zt = (Short_t)nucleus->getZ();
    }
    else {
      theEventInfo.At = (Short_t)nucleus->getA();
      theEventInfo.Zt = (Short_t)nucleus->getZ();
    }
    // Do nothing below the Coulomb barrier
    if(maxImpactParameter<=0.) {
      // Fill in the event information
    //Particle *pbar = new Particle;
    //PbarAtrestEntryChannel *obj = new PbarAtrestEntryChannel(nucleus, pbar);
      if(projectileSpecies.theType == antiProton && kineticEnergy <= theConfig->getAtrestThreshold()){         //D
        INCL_DEBUG("at rest annihilation" << '\n');
        //theEventInfo.transparent = false;
      } else {       
        theEventInfo.transparent = true;
        return false;
      }
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

  void INCL::postCascade(const ParticleSpecies &projectileSpecies, const G4double kineticEnergy) {
    // Fill in the event information
    theEventInfo.stoppingTime = propagationModel->getCurrentTime();

    // The event bias
    theEventInfo.eventBias = (Double_t) Particle::getTotalBias();

    // Forced CN?
    if(!(projectileSpecies.theType==antiProton && kineticEnergy<=theConfig->getAtrestThreshold())){
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
    }

    if(!(projectileSpecies.theType==antiProton && kineticEnergy<=theConfig->getAtrestThreshold())){
      theEventInfo.transparent = forceTransparent || nucleus->isEventTransparent();
    }

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
      theEventInfo.clusterDecay = nucleus->decayOutgoingClusters() || nucleus->decayMe(); //D

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
      theEventInfo.clusterDecay = nucleus->decayOutgoingClusters() || nucleus->decayMe(); //D

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
      else if(p.theType==antiProton) {
      maxUniverseRadius = rMax;                 //check interaction distance!!!
    }
    INCL_DEBUG("Initialised universe radius: " << maxUniverseRadius << '\n');
  }


  G4double INCL::initUniverseRadiusForAntiprotonAtRest(const G4int A, const G4int Z) {
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
    return rMax;                
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

  G4double INCL::read_file(std::string filename, std::vector<G4double>& probabilities, 
                           std::vector<std::vector<G4String>>& particle_types) {
    std::ifstream file(filename);
    G4double sum_probs = 0.0;
    if (file.is_open()) {
      G4String line;
      while (getline(file, line)) {
        std::istringstream iss(line);
        G4double prob;
        iss >> prob;
        sum_probs += prob;
        probabilities.push_back(prob);
        std::vector<G4String> types;
        G4String type;
        while (iss >> type) {
          types.push_back(type);
        }
        particle_types.push_back(types);
      }
    } else {
#ifdef INCLXX_IN_GEANT4_MODE
      G4cout << "ERROR no fread_file " << filename << G4endl;
#else
      std::cout << "ERROR no fread_file " << filename << std::endl;
#endif              
    }
    return sum_probs;
  }


  G4int INCL::findStringNumber(G4double rdm, std::vector<G4double> yields) {
    G4int stringNumber = -1;
    G4double smallestsum = 0.0;
    G4double biggestsum = yields[0];
    //G4cout << "initial input " << rdm << G4endl;
    for (G4int i = 0; i < static_cast<G4int>(yields.size()-1); i++) {
      if (rdm >= smallestsum && rdm <= biggestsum) {
        //G4cout << smallestsum << " and " << biggestsum << G4endl;
        stringNumber = i+1;
      }
      smallestsum += yields[i];
      biggestsum += yields[i+1];
    }
    if(stringNumber==-1) stringNumber = static_cast<G4int>(yields.size());
    if(stringNumber==-1){
      INCL_ERROR("ERROR in findStringNumber (stringNumber=-1)");
#ifdef INCLXX_IN_GEANT4_MODE
      G4cout << "ERROR in findStringNumber" << G4endl;
#else
      std::cout << "ERROR in findStringNumber" << std::endl;
#endif              
    }
    return stringNumber;
  }


  void INCL::preCascade_pbarH1(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy) {
    // Reset theEventInfo
    theEventInfo.reset();

    EventInfo::eventNumber++;

    // Fill in the event information
    theEventInfo.projectileType = projectileSpecies.theType;
    theEventInfo.Ap = -1;
    theEventInfo.Zp = -1;
    theEventInfo.Sp = 0;
    theEventInfo.Ep = kineticEnergy;
    theEventInfo.St = 0;
    theEventInfo.At = 1;
    theEventInfo.Zt = 1;
  }


  void INCL::postCascade_pbarH1(ParticleList const &outgoingParticles) {
    theEventInfo.nParticles = 0;

    // Reset the remnant counter
    theEventInfo.nRemnants = 0;
    theEventInfo.history.clear();

    for(ParticleIter i=outgoingParticles.begin(), e=outgoingParticles.end(); i!=e; ++i ) {
      theEventInfo.A[theEventInfo.nParticles] = (Short_t)(*i)->getA();
      theEventInfo.Z[theEventInfo.nParticles] = (Short_t)(*i)->getZ();
      theEventInfo.S[theEventInfo.nParticles] = (Short_t)(*i)->getS();
      theEventInfo.EKin[theEventInfo.nParticles] = (*i)->getKineticEnergy();
      ThreeVector mom = (*i)->getMomentum();
      theEventInfo.px[theEventInfo.nParticles] = mom.getX();
      theEventInfo.py[theEventInfo.nParticles] = mom.getY();
      theEventInfo.pz[theEventInfo.nParticles] = mom.getZ();
      theEventInfo.theta[theEventInfo.nParticles] = Math::toDegrees(mom.theta());
      theEventInfo.phi[theEventInfo.nParticles] = Math::toDegrees(mom.phi());
      theEventInfo.origin[theEventInfo.nParticles] = -1;
#ifdef INCLXX_IN_GEANT4_MODE
      theEventInfo.parentResonancePDGCode[theEventInfo.nParticles] = (*i)->getParentResonancePDGCode();
      theEventInfo.parentResonanceID[theEventInfo.nParticles] = (*i)->getParentResonanceID();
#endif
      theEventInfo.history.push_back("");
      ParticleSpecies pt((*i)->getType());
      theEventInfo.PDGCode[theEventInfo.nParticles] = pt.getPDGCode();
      theEventInfo.nParticles++;
    }
    theEventInfo.nCascadeParticles = theEventInfo.nParticles;
  }


  void INCL::preCascade_pbarH2(ParticleSpecies const &projectileSpecies, const G4double kineticEnergy) {
    // Reset theEventInfo
    theEventInfo.reset();

    EventInfo::eventNumber++;

    // Fill in the event information
    theEventInfo.projectileType = projectileSpecies.theType;
    theEventInfo.Ap = -1;
    theEventInfo.Zp = -1;
    theEventInfo.Sp = 0;
    theEventInfo.Ep = kineticEnergy;
    theEventInfo.St = 0;
    theEventInfo.At = 2;
    theEventInfo.Zt = 1;
  }


  void INCL::postCascade_pbarH2(ParticleList const &outgoingParticles, ParticleList const &H2Particles) {
    theEventInfo.nParticles = 0;

    // Reset the remnant counter
    theEventInfo.nRemnants = 0;
    theEventInfo.history.clear();

    for(ParticleIter i=outgoingParticles.begin(), e=outgoingParticles.end(); i!=e; ++i ) {
      theEventInfo.A[theEventInfo.nParticles] = (Short_t)(*i)->getA();
      theEventInfo.Z[theEventInfo.nParticles] = (Short_t)(*i)->getZ();
      theEventInfo.S[theEventInfo.nParticles] = (Short_t)(*i)->getS();
      theEventInfo.EKin[theEventInfo.nParticles] = (*i)->getKineticEnergy();
      ThreeVector mom = (*i)->getMomentum();
      theEventInfo.px[theEventInfo.nParticles] = mom.getX();
      theEventInfo.py[theEventInfo.nParticles] = mom.getY();
      theEventInfo.pz[theEventInfo.nParticles] = mom.getZ();
      theEventInfo.theta[theEventInfo.nParticles] = Math::toDegrees(mom.theta());
      theEventInfo.phi[theEventInfo.nParticles] = Math::toDegrees(mom.phi());
      theEventInfo.origin[theEventInfo.nParticles] = -1;
#ifdef INCLXX_IN_GEANT4_MODE
      theEventInfo.parentResonancePDGCode[theEventInfo.nParticles] = (*i)->getParentResonancePDGCode();
      theEventInfo.parentResonanceID[theEventInfo.nParticles] = (*i)->getParentResonanceID();
#endif
      theEventInfo.history.push_back("");
      ParticleSpecies pt((*i)->getType());
      theEventInfo.PDGCode[theEventInfo.nParticles] = pt.getPDGCode();
      theEventInfo.nParticles++;
    }

    for(ParticleIter i=H2Particles.begin(), e=H2Particles.end(); i!=e; ++i ) {
      theEventInfo.A[theEventInfo.nParticles] = (Short_t)(*i)->getA();
      theEventInfo.Z[theEventInfo.nParticles] = (Short_t)(*i)->getZ();
      theEventInfo.S[theEventInfo.nParticles] = (Short_t)(*i)->getS();
      theEventInfo.EKin[theEventInfo.nParticles] = (*i)->getKineticEnergy();
      ThreeVector mom = (*i)->getMomentum();
      theEventInfo.px[theEventInfo.nParticles] = mom.getX();
      theEventInfo.py[theEventInfo.nParticles] = mom.getY();
      theEventInfo.pz[theEventInfo.nParticles] = mom.getZ();
      theEventInfo.theta[theEventInfo.nParticles] = Math::toDegrees(mom.theta());
      theEventInfo.phi[theEventInfo.nParticles] = Math::toDegrees(mom.phi());
      theEventInfo.origin[theEventInfo.nParticles] = -1;
#ifdef INCLXX_IN_GEANT4_MODE
      theEventInfo.parentResonancePDGCode[theEventInfo.nParticles] = (*i)->getParentResonancePDGCode();
      theEventInfo.parentResonanceID[theEventInfo.nParticles] = (*i)->getParentResonanceID();
#endif
      theEventInfo.history.push_back("");
      ParticleSpecies pt((*i)->getType());
      theEventInfo.PDGCode[theEventInfo.nParticles] = pt.getPDGCode();
      theEventInfo.nParticles++;
    }
    theEventInfo.nCascadeParticles = theEventInfo.nParticles;
  }

}
