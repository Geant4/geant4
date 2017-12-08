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

#include "G4INCLParticleType.hh"
#include "G4INCLConfig.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {
  
  Config::Config() {
    init();
  }

  Config::~Config() {}

  void Config::init() {
    verbosity = 1;
    logFileName = "-";
    inputFileName = "";
    title = "INCL default run title";
    nShots = 1000;
    naturalTarget = false;
    projectileString = "proton";
    projectileSpecies = G4INCL::Proton;
    projectileKineticEnergy = 1000.0;
    verboseEvent = -1;
    randomSeeds = "";
    randomSeedVector.push_back(666);
    randomSeedVector.push_back(777);
    randomSeedVector.push_back(1234);
    pauliString = "strict-statistical";
    pauliType = StrictStatisticalPauli;
    CDPP = true;
    coulombString = "non-relativistic";
    coulombType = NonRelativisticCoulomb;
    potentialString = "isospin-energy";
    potentialType = IsospinEnergyPotential;
    pionPotential = true;
    localEnergyBBString = "first-collision";
    localEnergyBBType = FirstCollisionLocalEnergy;
    localEnergyPiString = "first-collision";
    localEnergyPiType = FirstCollisionLocalEnergy;
    deExcitationString = "none";
    deExcitationType = DeExcitationNone;
    clusterAlgorithmString = "intercomparison";
    clusterAlgorithmType = IntercomparisonClusterAlgorithm;
    clusterMaxMass = 8;
    backToSpectator = true;
    useRealMasses = true;
    impactParameter = -1.;
    separationEnergyString = "INCL";
    separationEnergyType = INCLSeparationEnergy;
    fermiMomentumString = "constant";
    fermiMomentumType = ConstantFermiMomentum;
    fermiMomentum = -1.;
    cutNN = 1910.;
#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
    maxMassFermiBreakUp = 16;
    maxChargeFermiBreakUp = 8;
#endif
    rpCorrelationCoefficient = 0.98;
    rpCorrelationCoefficientProton = 0.5;
    rpCorrelationCoefficientNeutron = 0.74;
    neutronSkin = 0.;
    neutronHalo = 0.;
    refraction=false;
    phaseSpaceGenerator = "Raubold-Lynch";
    phaseSpaceGeneratorType = RauboldLynchType;
    cascadeAction = "default";
    cascadeActionType = DefaultActionType;
    randomNumberGenerator = "Ranecu";
    rngType = RanecuType;
    autosaveFrequency = 10000;
    maxNumberMultipions = -1;
    crossSectionsString = "strangeness";
    crossSectionsType = StrangenessCrossSections;
    hadronizationTime = 0.;
#ifdef INCL_ROOT_USE
    conciseROOTTree = false;
#endif
    inverseKinematics = false;
    decayTimeThreshold = 1.e-20;
    bias = 1.;
  }

  std::string Config::summary() {
    std::stringstream message;
    message << "INCL++ version " << getVersionString() << '\n';
    if(projectileSpecies.theType != Composite)
      message << "Projectile: " << ParticleTable::getName(projectileSpecies) << '\n';
    else
      message << "Projectile: composite, A=" << projectileSpecies.theA << ", Z=" << projectileSpecies.theZ << ", S=" << projectileSpecies.theS << '\n';
    message << "  energy = " << projectileKineticEnergy << '\n';
    if(targetSpecies.theA>0)
      message << "Target: A = " << targetSpecies.theA << " Z = " << targetSpecies.theZ << " S = " << targetSpecies.theS << '\n';
    else
      message << "Target: natural isotopic composition, Z = " << targetSpecies.theZ << '\n';
    message << "Number of requested shots = " << nShots << '\n';
    return message.str();
  }

}
