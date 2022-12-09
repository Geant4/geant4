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

#include "G4INCLParticleTable.hh"
#include "G4INCLNuclearMassTable.hh"
#include <algorithm>
// #include <cassert>
#include <cmath>
#include <cctype>
#include <sstream>
#ifdef INCLXX_IN_GEANT4_MODE
#include "G4SystemOfUnits.hh"
#endif

#ifdef INCLXX_IN_GEANT4_MODE
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#endif

namespace G4INCL {

  namespace ParticleTable {

    namespace {

      /// \brief Static instance of the NaturalIsotopicAbundances class
      const NaturalIsotopicDistributions *theNaturalIsotopicDistributions = NULL;

      const G4double theINCLNucleonMass = 938.2796;
      const G4double theINCLPionMass = 138.0;
      const G4double theINCLLambdaMass = 1115.683;
//      const G4double theINCLSigmaMass = 1197.45;
//      const G4double theINCLKaonMass = 497.614;
      const G4double theINCLEtaMass = 547.862;
      const G4double theINCLOmegaMass = 782.65;
      const G4double theINCLEtaPrimeMass = 957.78;
      const G4double theINCLPhotonMass = 0.0;
      G4ThreadLocal G4double protonMass = 0.0;
      G4ThreadLocal G4double neutronMass = 0.0;
      G4ThreadLocal G4double piPlusMass = 0.0;
      G4ThreadLocal G4double piMinusMass = 0.0;
      G4ThreadLocal G4double piZeroMass = 0.0;
      G4ThreadLocal G4double SigmaPlusMass = 0.0;
      G4ThreadLocal G4double SigmaZeroMass = 0.0;
      G4ThreadLocal G4double SigmaMinusMass = 0.0;
      G4ThreadLocal G4double LambdaMass = 0.0;
      G4ThreadLocal G4double KPlusMass = 0.0;
      G4ThreadLocal G4double KZeroMass = 0.0;
      G4ThreadLocal G4double KZeroBarMass = 0.0;
      G4ThreadLocal G4double KShortMass = 0.0;
      G4ThreadLocal G4double KLongMass = 0.0;
      G4ThreadLocal G4double KMinusMass = 0.0;
      G4ThreadLocal G4double etaMass = 0.0;
      G4ThreadLocal G4double omegaMass = 0.0;
      G4ThreadLocal G4double etaPrimeMass = 0.0;
      G4ThreadLocal G4double photonMass = 0.0;

      // Hard-coded values of the real particle masses (MeV/c^2)
      G4ThreadLocal G4double theRealProtonMass = 938.27203;
      G4ThreadLocal G4double theRealNeutronMass = 939.56536;
      G4ThreadLocal G4double theRealChargedPiMass = 139.57018;
      G4ThreadLocal G4double theRealPiZeroMass = 134.9766;
      G4ThreadLocal G4double theRealLambdaMass = 1115.683;
      G4ThreadLocal G4double theRealSigmaPlusMass = 1189.37;
      G4ThreadLocal G4double theRealSigmaZeroMass = 1192.64;
      G4ThreadLocal G4double theRealSigmaMinusMass = 1197.45;
      G4ThreadLocal G4double theRealChargedKaonMass = 493.677;
      G4ThreadLocal G4double theRealNeutralKaonMass = 497.614;
      G4ThreadLocal G4double theRealEtaMass = 547.862;
      G4ThreadLocal G4double theRealOmegaMass = 782.65;
      G4ThreadLocal G4double theRealEtaPrimeMass = 957.78;
      G4ThreadLocal G4double theRealPhotonMass = 0.0;

      // Width (second)
      const G4double theChargedPiWidth = 2.6033e-08;
      const G4double thePiZeroWidth = 8.52e-17;
      const G4double theEtaWidth = 5.025e-19; // 1.31 keV
      const G4double theOmegaWidth = 7.7528e-23; // 8.49 MeV
      const G4double theEtaPrimeWidth = 3.3243e-21; // 0.198 MeV
      const G4double theChargedKaonWidth = 1.238e-08;
      const G4double theKShortWidth = 8.954e-11;
      const G4double theKLongWidth = 5.116e-08;
      const G4double theLambdaWidth = 2.632e-10;
      const G4double theSigmaPlusWidth = 8.018e-11;
      const G4double theSigmaZeroWidth = 7.4e-20;
      const G4double theSigmaMinusWidth = 1.479e-10;
      G4ThreadLocal G4double piPlusWidth = 0.0;
      G4ThreadLocal G4double piMinusWidth = 0.0;
      G4ThreadLocal G4double piZeroWidth = 0.0;
      G4ThreadLocal G4double etaWidth = 0.0;
      G4ThreadLocal G4double omegaWidth = 0.0;
      G4ThreadLocal G4double etaPrimeWidth = 0.0;
      G4ThreadLocal G4double LambdaWidth = 0.0;
      G4ThreadLocal G4double SigmaPlusWidth = 0.0;
      G4ThreadLocal G4double SigmaZeroWidth = 0.0;
      G4ThreadLocal G4double SigmaMinusWidth = 0.0;
      G4ThreadLocal G4double KPlusWidth = 0.0;
      G4ThreadLocal G4double KMinusWidth = 0.0;
      G4ThreadLocal G4double KShortWidth = 0.0;
      G4ThreadLocal G4double KLongWidth = 0.0;
        

      const G4int mediumNucleiTableSize = 30;

      const G4double mediumDiffuseness[mediumNucleiTableSize] =
      {0.0,0.0,0.0,0.0,0.0,1.78,1.77,1.77,1.69,1.71,
        1.69,1.72,1.635,1.730,1.81,1.833,1.798,
        1.93,0.567,0.571, 0.560,0.549,0.550,0.551,
        0.580,0.575,0.569,0.537,0.0,0.0};
      const G4double mediumRadius[mediumNucleiTableSize] =
      {0.0,0.0,0.0,0.0,0.0,0.334,0.327,0.479,0.631,0.838,
        0.811,0.84,1.403,1.335,1.25,1.544,1.498,1.57,
        2.58,2.77, 2.775,2.78,2.88,2.98,3.22,3.03,2.84,
        3.14,0.0,0.0};

      const G4double positionRMS[clusterTableZSize][clusterTableASize] = {
        /*     A=   0     1     2     3     4     5     6     7     8     9    10    11    12 */
        /* Z=0 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        /* Z=1 */ {-1.0, -1.0, 2.10, 1.80, 1.70, 1.83, 2.60, 2.50, -1.0, -1.0, -1.0, -1.0, -1.0},
        /* Z=2 */ {-1.0, -1.0, -1.0, 1.80, 1.68, 1.70, 2.60, 2.50, 2.50, 2.50, 2.50, -1.0, -1.0},
        /* Z=3 */ {-1.0, -1.0, -1.0, -1.0, 1.70, 1.83, 2.56, 2.40, 2.50, 2.50, 2.50, 2.50, 2.50},
        /* Z=4 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 2.60, 2.50, 2.50, 2.51, 2.50, 2.50, 2.50},
        /* Z=5 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 2.50, 2.50, 2.50, 2.50, 2.45, 2.40, 2.50},
        /* Z=6 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 2.50, 2.50, 2.50, 2.50, 2.47},
        /* Z=7 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 2.50, 2.50, 2.50},
        /* Z=8 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 2.50}
      };

      const G4double momentumRMS[clusterTableZSize][clusterTableASize] = {
        /*     A=   0     1     2     3     4     5     6     7     8     9    10    11    12 */
        /* Z=0 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0},
        /* Z=1 */ {-1.0, -1.0, 77.0, 110., 153., 100., 100., 100., -1.0, -1.0, -1.0, -1.0, -1.0},
        /* Z=2 */ {-1.0, -1.0, -1.0, 110., 153., 100., 100., 100., 100., 100., 100., -1.0, -1.0},
        /* Z=3 */ {-1.0, -1.0, -1.0, -1.0, 153., 100., 100., 100., 100., 100., 100., 100., 100.},
        /* Z=4 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 100., 100., 100., 100., 100., 100., 100.},
        /* Z=5 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 100., 100., 100., 100., 100., 100., 100.},
        /* Z=6 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 100., 100., 100., 100., 100.},
        /* Z=7 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 100., 100., 100.},
        /* Z=8 */ {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 100.}
      };

      const G4int elementTableSize = 113; // up to Cn

      /// \brief Table of chemical element names
      const std::string elementTable[elementTableSize] = {
        "",
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
        "Cn"
      };

      /// \brief Digit names to compose IUPAC element names
      const std::string elementIUPACDigits = "nubtqphsoe";

#define INCL_DEFAULT_SEPARATION_ENERGY 6.83
      const G4double theINCLProtonSeparationEnergy = INCL_DEFAULT_SEPARATION_ENERGY;
      const G4double theINCLNeutronSeparationEnergy = INCL_DEFAULT_SEPARATION_ENERGY;
      const G4double theINCLLambdaSeparationEnergy = INCL_DEFAULT_SEPARATION_ENERGY;
      G4ThreadLocal G4double protonSeparationEnergy = INCL_DEFAULT_SEPARATION_ENERGY;
      G4ThreadLocal G4double neutronSeparationEnergy = INCL_DEFAULT_SEPARATION_ENERGY;
      G4ThreadLocal G4double lambdaSeparationEnergy = INCL_DEFAULT_SEPARATION_ENERGY;
#undef INCL_DEFAULT_SEPARATION_ENERGY

      G4ThreadLocal G4double rpCorrelationCoefficient[UnknownParticle];

      G4ThreadLocal G4double neutronSkin = 0.0;
      G4ThreadLocal G4double neutronHalo = 0.0;

#ifdef INCLXX_IN_GEANT4_MODE
      G4ThreadLocal G4IonTable *theG4IonTable;
#endif

      /// \brief Default value for constant Fermi momentum
      G4ThreadLocal G4double constantFermiMomentum = 0.0;

      /// \brief Transform a IUPAC char to an char representing an integer digit
      char iupacToInt(char c) {
        return (char)(((G4int)'0')+elementIUPACDigits.find(c));
      }

      /// \brief Transform an integer digit (represented by a char) to a IUPAC char
      char intToIUPAC(char n) { return elementIUPACDigits.at(n); }

      /// \brief Get the singleton instance of the natural isotopic distributions
      const NaturalIsotopicDistributions *getNaturalIsotopicDistributions() {
        if(!theNaturalIsotopicDistributions)
          theNaturalIsotopicDistributions = new NaturalIsotopicDistributions;
        return theNaturalIsotopicDistributions;
      }

    } // namespace

    void initialize(Config const * const theConfig /*=0*/) {
      protonMass = theINCLNucleonMass;
      neutronMass = theINCLNucleonMass;
      piPlusMass = theINCLPionMass;
      piMinusMass = theINCLPionMass;
      piZeroMass = theINCLPionMass;
      /*
      SigmaPlusMass = theINCLSigmaMass;
      SigmaMinusMass = theINCLSigmaMass;
      SigmaZeroMass = theINCLSigmaMass;
      LambdaMass = theINCLLambdaMass;
      KPlusMass = theINCLKaonMass;
      KZeroMass = theINCLKaonMass;
      KZeroBarMass = theINCLKaonMass;
      KShortMass = theINCLKaonMass;
      KLongMass = theINCLKaonMass;
      KMinusMass = theINCLKaonMass;
      */
      SigmaPlusMass = theRealSigmaPlusMass;
      SigmaMinusMass = theRealSigmaMinusMass;
      SigmaZeroMass = theRealSigmaZeroMass;
      LambdaMass = theINCLLambdaMass;
      KPlusMass = theRealChargedKaonMass;
      KZeroMass = theRealNeutralKaonMass;
      KZeroBarMass = theRealNeutralKaonMass;
      KShortMass = theRealNeutralKaonMass;
      KLongMass = theRealNeutralKaonMass;
      KMinusMass = theRealChargedKaonMass;
      
      etaMass = theINCLEtaMass;
      omegaMass = theINCLOmegaMass;
      etaPrimeMass = theINCLEtaPrimeMass;
      photonMass = theINCLPhotonMass;
      if(theConfig && theConfig->getUseRealMasses()) {
        getTableMass = getRealMass;
        getTableParticleMass = getRealMass;
      } else {
        getTableMass = getINCLMass;
        getTableParticleMass = getINCLMass;
      }

#ifndef INCLXX_IN_GEANT4_MODE
      std::string dataFilePath;
      if(theConfig)
        dataFilePath = theConfig->getINCLXXDataFilePath();
      NuclearMassTable::initialize(dataFilePath, getRealMass(Proton), getRealMass(Neutron));
#endif

#ifdef INCLXX_IN_GEANT4_MODE
      G4ParticleTable *theG4ParticleTable = G4ParticleTable::GetParticleTable();
      theG4IonTable = theG4ParticleTable->GetIonTable();
      theRealProtonMass = theG4ParticleTable->FindParticle("proton")->GetPDGMass() / MeV;
      theRealNeutronMass = theG4ParticleTable->FindParticle("neutron")->GetPDGMass() / MeV;
      theRealChargedPiMass = theG4ParticleTable->FindParticle("pi+")->GetPDGMass() / MeV;
      theRealPiZeroMass = theG4ParticleTable->FindParticle("pi0")->GetPDGMass() / MeV;
      theRealEtaMass = theG4ParticleTable->FindParticle("eta")->GetPDGMass() / MeV;
      theRealOmegaMass = theG4ParticleTable->FindParticle("omega")->GetPDGMass() / MeV;
      theRealEtaPrimeMass = theG4ParticleTable->FindParticle("eta_prime")->GetPDGMass() / MeV;
      theRealPhotonMass = theG4ParticleTable->FindParticle("gamma")->GetPDGMass() / MeV;
      theRealSigmaPlusMass = theG4ParticleTable->FindParticle("sigma+")->GetPDGMass() / MeV;
      theRealSigmaZeroMass = theG4ParticleTable->FindParticle("sigma0")->GetPDGMass() / MeV;
      theRealSigmaMinusMass = theG4ParticleTable->FindParticle("sigma-")->GetPDGMass() / MeV;
      theRealLambdaMass = theG4ParticleTable->FindParticle("lambda")->GetPDGMass() / MeV;
      theRealChargedKaonMass = theG4ParticleTable->FindParticle("kaon+")->GetPDGMass() / MeV;
      theRealNeutralKaonMass = theG4ParticleTable->FindParticle("kaon0")->GetPDGMass() / MeV;
#endif

      minDeltaMass = theRealNeutronMass + theRealChargedPiMass + 0.5;
      minDeltaMass2 = minDeltaMass*minDeltaMass;
      minDeltaMassRndm = std::atan((minDeltaMass-effectiveDeltaMass)*2./effectiveDeltaWidth);

      piPlusWidth   = theChargedPiWidth;
      piMinusWidth  = theChargedPiWidth;
      piZeroWidth   = thePiZeroWidth;
      etaWidth      = theEtaWidth;
      omegaWidth    = theOmegaWidth;
      etaPrimeWidth = theEtaPrimeWidth;
      
      SigmaMinusWidth = theSigmaMinusWidth;
      SigmaPlusWidth = theSigmaPlusWidth;
      SigmaZeroWidth = theSigmaZeroWidth;
      LambdaWidth = theLambdaWidth;
      KPlusWidth = theChargedKaonWidth;
      KMinusWidth = theChargedKaonWidth;
      KShortWidth = theKShortWidth;
      KLongWidth = theKLongWidth;
        
      // Initialise HFB tables
#ifdef INCLXX_IN_GEANT4_MODE
        HFB::initialize();
#else
        HFB::initialize(dataFilePath);
#endif

      // Initialise the separation-energy function
      if(!theConfig || theConfig->getSeparationEnergyType()==INCLSeparationEnergy)
        getSeparationEnergy = getSeparationEnergyINCL;
      else if(theConfig->getSeparationEnergyType()==RealSeparationEnergy)
        getSeparationEnergy = getSeparationEnergyReal;
      else if(theConfig->getSeparationEnergyType()==RealForLightSeparationEnergy)
        getSeparationEnergy = getSeparationEnergyRealForLight;
      else {
        INCL_FATAL("Unrecognized separation-energy type in ParticleTable initialization: " << theConfig->getSeparationEnergyType() << '\n');
        return;
      }

      // Initialise the Fermi-momentum function
      if(!theConfig || theConfig->getFermiMomentumType()==ConstantFermiMomentum) {
        getFermiMomentum = ParticleTable::getFermiMomentumConstant;
        if(theConfig) {
          const G4double aFermiMomentum = theConfig->getFermiMomentum();
          if(aFermiMomentum>0.)
            constantFermiMomentum = aFermiMomentum;
          else
            constantFermiMomentum = PhysicalConstants::Pf;
        } else {
          constantFermiMomentum = PhysicalConstants::Pf;
        }
      } else if(theConfig->getFermiMomentumType()==ConstantLightFermiMomentum)
        getFermiMomentum = ParticleTable::getFermiMomentumConstantLight;
      else if(theConfig->getFermiMomentumType()==MassDependentFermiMomentum)
        getFermiMomentum = ParticleTable::getFermiMomentumMassDependent;
      else {
        INCL_FATAL("Unrecognized Fermi-momentum type in ParticleTable initialization: " << theConfig->getFermiMomentumType() << '\n');
        return;
      }

      // Initialise the r-p correlation coefficients
      std::fill(rpCorrelationCoefficient, rpCorrelationCoefficient + UnknownParticle, 1.);
      if(theConfig) {
        rpCorrelationCoefficient[Proton] = theConfig->getRPCorrelationCoefficient(Proton);
        rpCorrelationCoefficient[Neutron] = theConfig->getRPCorrelationCoefficient(Neutron);
      }

      // Initialise the neutron-skin parameters
      if(theConfig) {
        neutronSkin = theConfig->getNeutronSkin();
        neutronHalo = theConfig->getNeutronHalo();
      }

    }

    G4int getIsospin(const ParticleType t) {
      // Actually this is the 3rd component of isospin (I_z) multiplied by 2!
      if(t == Proton) {
        return 1;
      } else if(t == Neutron) {
        return -1;
      } else if(t == PiPlus) {
        return 2;
      } else if(t == PiMinus) {
        return -2;
      } else if(t == PiZero) {
        return 0;
      } else if(t == DeltaPlusPlus) {
        return 3;
      } else if(t == DeltaPlus) {
        return 1;
      } else if(t == DeltaZero) {
        return -1;
      } else if(t == DeltaMinus) {
        return -3;
      } else if(t == Lambda) {
        return 0;
      } else if(t == SigmaPlus) {
        return 2;
      } else if(t == SigmaZero) {
        return 0;
      } else if(t == SigmaMinus) {
        return -2;
      } else if(t == KPlus) {
        return 1;
      } else if(t == KZero) {
        return -1;
      } else if(t == KZeroBar) {
        return 1;
      } else if(t == KShort) {
        return 0;
      } else if(t == KLong) {
        return 0;
      } else if(t == KMinus) {
        return -1;
      } else if(t == Eta) {
        return 0;
      } else if(t == Omega) {
        return 0;
      } else if(t == EtaPrime) {
        return 0;
      } else if(t == Photon) {
        return 0;
      }
      INCL_ERROR("Requested isospin of an unknown particle!");
      return -10; // Unknown
    }

    std::string getShortName(const ParticleSpecies &sp) {
      if(sp.theType==Composite && sp.theS == 0)
        return getShortName(sp.theA,sp.theZ);
      else if(sp.theType==Composite)
        return getName(sp.theA,sp.theZ,sp.theS);
      else
        return getShortName(sp.theType);
    }

    std::string getName(const ParticleSpecies &sp) {
      if(sp.theType==Composite && sp.theS == 0)
        return getName(sp.theA,sp.theZ);
      else if(sp.theType==Composite)
        return getName(sp.theA,sp.theZ,sp.theS);
      else
        return getName(sp.theType);
    }

    std::string getName(const G4int A, const G4int Z) {
      std::stringstream stream;
      stream << getElementName(Z) << "-" << A;
      return stream.str();
    }

    std::string getName(const G4int A, const G4int Z, const G4int S) {
      std::stringstream stream;
      if(S >= 0) // S < 0 for hypernuclei
        return getName(A, Z);
      else if(S == -1)
        stream << getElementName(Z) << "-" << A << "_" << "Lambda";
      else
        stream << getElementName(Z) << "-" << A << "_" << S << "-Lambda";
      return stream.str();
    }

    std::string getShortName(const G4int A, const G4int Z) {
      std::stringstream stream;
      stream << getElementName(Z);
      if(A>0)
        stream << A;
      return stream.str();
    }

    std::string getName(const ParticleType p) {
      if(p == G4INCL::Proton) {
        return std::string("proton");
      } else if(p == G4INCL::Neutron) {
        return std::string("neutron");
      } else if(p == G4INCL::DeltaPlusPlus) {
        return std::string("delta++");
      } else if(p == G4INCL::DeltaPlus) {
        return std::string("delta+");
      } else if(p == G4INCL::DeltaZero) {
        return std::string("delta0");
      } else if(p == G4INCL::DeltaMinus) {
        return std::string("delta-");
      } else if(p == G4INCL::PiPlus) {
        return std::string("pi+");
      } else if(p == G4INCL::PiZero) {
        return std::string("pi0");
      } else if(p == G4INCL::PiMinus) {
        return std::string("pi-");
      } else if(p == G4INCL::Lambda) {
        return std::string("lambda");
      } else if(p == G4INCL::SigmaPlus) {
        return std::string("sigma+");
      } else if(p == G4INCL::SigmaZero) {
        return std::string("sigma0");
      } else if(p == G4INCL::SigmaMinus) {
        return std::string("sigma-");
      } else if(p == G4INCL::KPlus) {
        return std::string("kaon+");
      } else if(p == G4INCL::KZero) {
        return std::string("kaon0");
      } else if(p == G4INCL::KZeroBar) {
        return std::string("kaon0bar");
      } else if(p == G4INCL::KMinus) {
        return std::string("kaon-");
      } else if(p == G4INCL::KShort) {
        return std::string("kaonshort");
      } else if(p == G4INCL::KLong) {
        return std::string("kaonlong");
      } else if(p == G4INCL::Composite) {
        return std::string("composite");
      } else if(p == G4INCL::Eta) {
        return std::string("eta");
      } else if(p == G4INCL::Omega) {
        return std::string("omega");
      } else if(p == G4INCL::EtaPrime) {
        return std::string("etaprime");
      } else if(p == G4INCL::Photon) {
        return std::string("photon");
      }
      return std::string("unknown");
    }

    std::string getShortName(const ParticleType p) {
      if(p == G4INCL::Proton) {
        return std::string("p");
      } else if(p == G4INCL::Neutron) {
        return std::string("n");
      } else if(p == G4INCL::DeltaPlusPlus) {
        return std::string("d++");
      } else if(p == G4INCL::DeltaPlus) {
        return std::string("d+");
      } else if(p == G4INCL::DeltaZero) {
        return std::string("d0");
      } else if(p == G4INCL::DeltaMinus) {
        return std::string("d-");
      } else if(p == G4INCL::PiPlus) {
        return std::string("pi+");
      } else if(p == G4INCL::PiZero) {
        return std::string("pi0");
      } else if(p == G4INCL::PiMinus) {
        return std::string("pi-");
      } else if(p == G4INCL::Lambda) {
        return std::string("l");
      } else if(p == G4INCL::SigmaPlus) {
        return std::string("s+");
      } else if(p == G4INCL::SigmaZero) {
        return std::string("s0");
      } else if(p == G4INCL::SigmaMinus) {
        return std::string("s-");
      } else if(p == G4INCL::KPlus) {
        return std::string("k+");
      } else if(p == G4INCL::KZero) {
        return std::string("k0");
      } else if(p == G4INCL::KZeroBar) {
        return std::string("k0b");
      } else if(p == G4INCL::KMinus) {
        return std::string("k-");
      } else if(p == G4INCL::KShort) {
        return std::string("ks");
      } else if(p == G4INCL::KLong) {
        return std::string("kl");
      } else if(p == G4INCL::Composite) {
        return std::string("comp");
      } else if(p == G4INCL::Eta) {
        return std::string("eta");
      } else if(p == G4INCL::Omega) {
        return std::string("omega");
      } else if(p == G4INCL::EtaPrime) {
        return std::string("etap");
      } else if(p == G4INCL::Photon) {
        return std::string("photon");
      }
      return std::string("unknown");
    }

    G4double getINCLMass(const ParticleType pt) {
      if(pt == Proton) {
        return protonMass;
      } else if(pt == Neutron) {
        return neutronMass;
      } else if(pt == PiPlus) {
        return piPlusMass;
      } else if(pt == PiMinus) {
        return piMinusMass;
      } else if(pt == PiZero) {
        return piZeroMass;
      } else if(pt == SigmaPlus) {
        return SigmaPlusMass;
      } else if(pt == SigmaMinus) {
        return SigmaMinusMass;
      } else if(pt == SigmaZero) {
        return SigmaZeroMass;
      } else if(pt == Lambda) {
        return LambdaMass;
      } else if(pt == KPlus) {
        return KPlusMass;
      } else if(pt == KZero) {
        return KZeroMass;
      } else if(pt == KZeroBar) {
        return KZeroBarMass;
      } else if(pt == KMinus) {
        return KMinusMass;
      } else if(pt == KShort) {
        return KShortMass;
      } else if(pt == KLong) {
        return KLongMass;
      } else if(pt == Eta) {
        return etaMass;
      } else if(pt == Omega) {
        return omegaMass;
      } else if(pt == EtaPrime) {
        return etaPrimeMass;
      } else if(pt == Photon) {
        return photonMass;
      } else {
        INCL_ERROR("getMass : Unknown particle type." << '\n');
        return 0.0;
      }
    }
          
    G4double getRealMass(const ParticleType t) {
      switch(t) {
        case Proton:
          return theRealProtonMass;
          break;
        case Neutron:
          return theRealNeutronMass;
          break;
        case PiPlus:
        case PiMinus:
          return theRealChargedPiMass;
          break;
        case PiZero:
          return theRealPiZeroMass;
          break;
        case SigmaPlus:
          return theRealSigmaPlusMass;
          break;
        case SigmaZero:
          return theRealSigmaZeroMass;
          break;
        case SigmaMinus:
          return theRealSigmaMinusMass;
          break;
        case Lambda:
          return theRealLambdaMass;
          break;
        case KPlus:
        case KMinus:
          return theRealChargedKaonMass;
          break;
        case KZero:
        case KZeroBar:
        case KShort:
        case KLong:
          return theRealNeutralKaonMass;
          break;
        case Eta:
          return theRealEtaMass;
          break;
        case Omega:
          return theRealOmegaMass;
          break;
        case EtaPrime:
          return theRealEtaPrimeMass;
          break;
        case Photon:
          return theRealPhotonMass;
          break;
        default:
          INCL_ERROR("Particle::getRealMass : Unknown particle type." << '\n');
          return 0.0;
          break;
      }
    }
    
    G4double getRealMass(const G4int A, const G4int Z, const G4int S) {
// assert(A>=0);
      // For nuclei with Z<0 or Z>A, assume that the exotic charge state is due to pions
      if(Z<0 && S<0)
        return (A+S)*theRealNeutronMass - S*LambdaMass - Z*getRealMass(PiMinus);
      else if(Z>A && S<0)
        return (A+S)*theRealProtonMass - S*LambdaMass + (A+S-Z)*getRealMass(PiPlus);
      if(Z<0)
        return (A)*theRealNeutronMass - Z*getRealMass(PiMinus);
      else if(Z>A)
        return (A)*theRealProtonMass + (A-Z)*getRealMass(PiPlus);
      else if(Z==0 && S==0)
        return A*theRealNeutronMass;
      else if(A==Z)
        return A*theRealProtonMass;
      else if(Z==0 && S<0)
        return (A+S)*theRealNeutronMass-S*LambdaMass;
      else if(A>1) {
#ifndef INCLXX_IN_GEANT4_MODE
        return ::G4INCL::NuclearMassTable::getMass(A,Z,S);
#else
        if(S<0) return theG4IonTable->GetNucleusMass(Z,A,std::abs(S)) / MeV;
        else    return theG4IonTable->GetNucleusMass(Z,A) / MeV;
#endif
      } else
        return 0.;
    }

    G4double getINCLMass(const G4int A, const G4int Z, const G4int S) {
// assert(A>=0);
      // For nuclei with Z<0 or Z>A, assume that the exotic charge state is due to pions
      // Note that S<0 for lambda
      if(Z<0 && S<0)
        return (A+S)*neutronMass - S*LambdaMass - Z*getINCLMass(PiMinus);
      else if(Z>A && S<0)
        return (A+S)*protonMass - S*LambdaMass + (A+S-Z)*getINCLMass(PiPlus);  
      else if(Z<0)
        return (A)*neutronMass - Z*getINCLMass(PiMinus);
      else if(Z>A)
        return (A)*protonMass + (A-Z)*getINCLMass(PiPlus);
      else if(A>1 && S<0)
        return Z*(protonMass - protonSeparationEnergy) + (A+S-Z)*(neutronMass - neutronSeparationEnergy) + std::abs(S)*(LambdaMass - lambdaSeparationEnergy);
      else if(A>1)
        return Z*(protonMass - protonSeparationEnergy) + (A-Z)*(neutronMass - neutronSeparationEnergy);
      else if(A==1 && Z==0 && S==0)
        return getINCLMass(Neutron);
      else if(A==1 && Z==1 && S==0)
        return getINCLMass(Proton);
      else if(A==1 && Z==0 && S==-1)
        return getINCLMass(Lambda);
      else
        return 0.;
    }

    G4double getTableQValue(const G4int A1, const G4int Z1, const G4int S1, const G4int A2, const G4int Z2, const G4int S2) {
      return getTableMass(A1,Z1,S1) + getTableMass(A2,Z2,S2) - getTableMass(A1+A2,Z1+Z2,S1+S2);
    }

    G4double getTableQValue(const G4int A1, const G4int Z1, const G4int S1, const G4int A2, const G4int Z2, const G4int S2, const G4int A3, const G4int Z3, const G4int S3) {
      return getTableMass(A1,Z1,S1) + getTableMass(A2,Z2,S2) - getTableMass(A3,Z3,S3) - getTableMass(A1+A2-A3,Z1+Z2-Z3,S1+S2-S3);
    }

    G4double getTableSpeciesMass(const ParticleSpecies &p) {
      if(p.theType == Composite)
        return (*getTableMass)(p.theA, p.theZ, p.theS);
      else
        return (*getTableParticleMass)(p.theType);
    }

    G4int getMassNumber(const ParticleType t) {
      switch(t) {
        case Proton:
        case Neutron:
        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
        case SigmaPlus:
        case SigmaZero:
        case SigmaMinus:
        case Lambda:
          return 1;
          break;
        case PiPlus:
        case PiMinus:
        case PiZero:
        case KPlus:
        case KZero:
        case KZeroBar:
        case KShort:
        case KLong:
        case KMinus:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          return 0;
          break;
        default:
          return 0;
          break;
      }
    }

    G4int getChargeNumber(const ParticleType t) {
      switch(t) {
        case DeltaPlusPlus:
          return 2;
          break;
        case Proton:
        case DeltaPlus:
        case PiPlus:
        case SigmaPlus:
        case KPlus:
          return 1;
          break;
        case Neutron:
        case DeltaZero:
        case PiZero:
        case SigmaZero:
        case Lambda:
        case KZero:
        case KZeroBar:
        case KShort:
        case KLong:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          return 0;
          break;
        case DeltaMinus:
        case PiMinus:
        case SigmaMinus:
        case KMinus:
          return -1;
          break;
        default:
          return 0;
          break;
      }
    }
    
    G4int getStrangenessNumber(const ParticleType t) {
      switch(t) {
        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
        case Proton:
        case Neutron:
        case PiPlus:
        case PiZero:
        case PiMinus:
        case Eta:
        case Omega:
        case EtaPrime:
        case Photon:
          return 0;
          break;
        case Lambda:
        case SigmaPlus:
        case SigmaZero:
        case SigmaMinus:
        case KZeroBar:
        case KMinus:
          return -1;
          break;
        case KPlus:
        case KZero:
          return 1;
          break;
        case KShort:
          return 0;
          break;
        case KLong:
          return 0;
          break;
        default:
          return 0;
          break;
      }
    }

    G4double getNuclearRadius(const ParticleType t, const G4int A, const G4int Z) {
// assert(A>=0);
      if(A > 19 || (A < 6 && A >= 2)) {
        // For large (Woods-Saxon or Modified Harmonic Oscillator) or small
        // (Gaussian) nuclei, the radius parameter is just the nuclear radius
        return getRadiusParameter(t,A,Z);
      } else if(A < clusterTableASize && Z>=0 && Z < clusterTableZSize && A >= 6) {
        const G4double thisRMS = positionRMS[Z][A];
        if(thisRMS>0.0)
          return thisRMS;
        else {
          INCL_DEBUG("getNuclearRadius: Radius for nucleus A = " << A << " Z = " << Z << " is not available" << '\n'
                     << "returning radius for C12");
          return positionRMS[6][12];
        }
      } else if(A <= 19) {
        const G4double theRadiusParameter = getRadiusParameter(t, A, Z);
        const G4double theDiffusenessParameter = getSurfaceDiffuseness(t, A, Z);
        // The formula yields the nuclear RMS radius based on the parameters of
        // the nuclear-density function
        return 1.225*theDiffusenessParameter*
          std::sqrt((2.+5.*theRadiusParameter)/(2.+3.*theRadiusParameter));
      } else {
        INCL_ERROR("getNuclearRadius: No radius for nucleus A = " << A << " Z = " << Z << '\n');
        return 0.0;
      }
    }

    G4double getLargestNuclearRadius(const G4int A, const G4int Z) {
      return Math::max(getNuclearRadius(Proton, A, Z), getNuclearRadius(Neutron, A, Z));
    }

    G4double getRadiusParameter(const ParticleType t, const G4int A, const G4int Z) {
// assert(A>0);
      if(A > 19) {
        // radius fit for lambdas
        if(t==Lambda){
         G4double r0 = (1.128+0.439*std::pow(A,-2./3.)) * std::pow(A, 1.0/3.0);
         return r0;
        }
        // phenomenological radius fit
        G4double r0 = (2.745e-4 * A + 1.063) * std::pow(A, 1.0/3.0);
        // HFB calculations
        if(getRPCorrelationCoefficient(t)<1.){
         G4double r0hfb = HFB::getRadiusParameterHFB(t,A,Z);
         if(r0hfb>0.)r0 = r0hfb;
        }
        //
        if(t==Neutron)
          r0 += neutronSkin;
        return r0;
      } else if(A < 6 && A >= 2) {
        if(Z<clusterTableZSize && Z>=0) {
          const G4double thisRMS = positionRMS[Z][A];
          if(thisRMS>0.0)
            return thisRMS;
          else {
            INCL_DEBUG("getRadiusParameter: Radius for nucleus A = " << A << " Z = " << Z << " is not available" << '\n'
                       << "returning radius for C12");
            return positionRMS[6][12];
          }
        } else {
          INCL_DEBUG("getRadiusParameter: Radius for nucleus A = " << A << " Z = " << Z << " is not available" << '\n'
                     << "returning radius for C12");
          return positionRMS[6][12];
        }
      } else if(A <= 19 && A >= 6) {
        if(t==Lambda){
         G4double r0 = (1.128+0.439*std::pow(A,-2./3.)) * std::pow(A, 1.0/3.0);
         return r0;
        }
        // HFB calculations
        if(getRPCorrelationCoefficient(t)<1.){
         G4double r0hfb = HFB::getSurfaceDiffusenessHFB(t,A,Z);
         if(r0hfb>0.)return r0hfb;
        }
        return mediumRadius[A-1];
        //      return 1.581*mediumDiffuseness[A-1]*(2.+5.*mediumRadius[A-1])/(2.+3.*mediumRadius[A-1]);
      } else {
        INCL_ERROR("getRadiusParameter: No radius for nucleus A = " << A << " Z = " << Z << '\n');
        return 0.0;
      }
    }

    G4double getMaximumNuclearRadius(const ParticleType t, const G4int A, const G4int Z) {
      const G4double XFOISA = 8.0;
      if(A > 19) {
        return getNuclearRadius(t,A,Z) + XFOISA * getSurfaceDiffuseness(t,A,Z);
      } else if(A <= 19 && A >= 6) {
        return 5.5 + 0.3 * (G4double(A) - 6.0)/12.0;
      } else if(A >= 2) {
        return getNuclearRadius(t, A, Z) + 4.5;
      } else {
        INCL_ERROR("getMaximumNuclearRadius : No maximum radius for nucleus A = " << A << " Z = " << Z << '\n');
        return 0.0;
      }
    }

    G4double getSurfaceDiffuseness(const ParticleType t, const G4int A, const G4int Z) {
      if(A > 19) {
        // phenomenological fit
        G4double a = 1.63e-4 * A + 0.510;
        // HFB calculations
        if(getRPCorrelationCoefficient(t)<1.){
          G4double ahfb = HFB::getSurfaceDiffusenessHFB(t,A,Z);
          if(ahfb>0.)a=ahfb;
        }
        //
        if(t==Lambda){
        // Like for neutrons
          G4double ahfb = HFB::getSurfaceDiffusenessHFB(Neutron,A,Z);
          if(ahfb>0.)a=ahfb;
        }
        if(t==Neutron)
          a += neutronHalo;
        return a;
      } else if(A <= 19 && A >= 6) {
        // HFB calculations
        if(getRPCorrelationCoefficient(t)<1.){
          G4double ahfb = HFB::getRadiusParameterHFB(t,A,Z);
          if(ahfb>0.)return ahfb;
        }
        return mediumDiffuseness[A-1];
      } else if(A < 6 && A >= 2) {
        INCL_ERROR("getSurfaceDiffuseness: was called for A = " << A << " Z = " << Z << '\n');
        return 0.0;
      } else {
        INCL_ERROR("getSurfaceDiffuseness: No diffuseness for nucleus A = " << A << " Z = " << Z << '\n');
        return 0.0;
      }
    }

    G4double getMomentumRMS(const G4int A, const G4int Z) {
// assert(Z>=0 && A>=0 && Z<=A);
      return getFermiMomentum(A,Z) * Math::sqrtThreeFifths;
    }

    G4double getSeparationEnergyINCL(const ParticleType t, const G4int /*A*/, const G4int /*Z*/) {
      if(t==Proton)
        return theINCLProtonSeparationEnergy;
      else if(t==Neutron)
        return theINCLNeutronSeparationEnergy;
      else if(t==Lambda)
        return theINCLLambdaSeparationEnergy;
      else {
        INCL_ERROR("ParticleTable::getSeparationEnergyINCL : Unknown particle type." << '\n');
        return 0.0;
      }
    }

    G4double getSeparationEnergyReal(const ParticleType t, const G4int A, const G4int Z) {
      // Real separation energies for all nuclei
      if(t==Proton)
        return (*getTableParticleMass)(Proton) + (*getTableMass)(A-1,Z-1,0) - (*getTableMass)(A,Z,0);
      else if(t==Neutron)
        return (*getTableParticleMass)(Neutron) + (*getTableMass)(A-1,Z,0) - (*getTableMass)(A,Z,0);
      else if(t==Lambda)
        return (*getTableParticleMass)(Lambda) + (*getTableMass)(A-1,Z,0) - (*getTableMass)(A,Z,-1);
      else {
        INCL_ERROR("ParticleTable::getSeparationEnergyReal : Unknown particle type." << '\n');
        return 0.0;
      }
    }

    G4double getSeparationEnergyRealForLight(const ParticleType t, const G4int A, const G4int Z) {
      // Real separation energies for light nuclei, fixed values for heavy nuclei
      if(Z<clusterTableZSize && A<clusterTableASize)
        return getSeparationEnergyReal(t, A, Z);
      else
        return getSeparationEnergyINCL(t, A, Z);
    }

    G4double getProtonSeparationEnergy() { return protonSeparationEnergy; }

    G4double getNeutronSeparationEnergy() { return neutronSeparationEnergy; }

    G4double getLambdaSeparationEnergy() { return lambdaSeparationEnergy; }

    void setProtonSeparationEnergy(const G4double sen) { protonSeparationEnergy = sen; }

    void setNeutronSeparationEnergy(const G4double sen) { neutronSeparationEnergy  = sen; }

    void setLambdaSeparationEnergy(const G4double sen) { lambdaSeparationEnergy  = sen; }

    std::string getElementName(const G4int Z) {
      if(Z<1) {
        INCL_WARN("getElementName called with Z<1" << '\n');
        return elementTable[0];
      } else if(Z<elementTableSize)
        return elementTable[Z];
      else
        return getIUPACElementName(Z);
    }

    std::string getIUPACElementName(const G4int Z) {
      std::stringstream elementStream;
      elementStream << Z;
      std::string elementName = elementStream.str();
      std::transform(elementName.begin(), elementName.end(), elementName.begin(), intToIUPAC);
      elementName[0] = (char)std::toupper(elementName.at(0));
      return elementName;
    }

    G4int parseElement(std::string pS) {
      // Normalize the element name
      std::transform(pS.begin(), pS.end(), pS.begin(), ::tolower);
      pS[0] = (char)std::toupper(pS[0]);

      const std::string *iter = std::find(elementTable, elementTable+elementTableSize, pS);
      if(iter != elementTable+elementTableSize)
        return G4int(iter - elementTable);
      else
        return ParticleTable::parseIUPACElement(pS);
    }

    G4int parseIUPACElement(std::string const &sel) {
      // Normalise to lower case
      std::string elementName(sel);
      std::transform(elementName.begin(), elementName.end(), elementName.begin(), ::tolower);
      // Return 0 if the element name contains anything but IUPAC digits
      if(elementName.find_first_not_of(elementIUPACDigits)!=std::string::npos)
        return 0;
      std::transform(elementName.begin(), elementName.end(), elementName.begin(), iupacToInt);
      std::stringstream elementStream(elementName);
      G4int Z;
      elementStream >> Z;
      return Z;
    }

    IsotopicDistribution const &getNaturalIsotopicDistribution(const G4int Z) {
      return getNaturalIsotopicDistributions()->getIsotopicDistribution(Z);
    }

    G4int drawRandomNaturalIsotope(const G4int Z) {
      return getNaturalIsotopicDistributions()->drawRandomIsotope(Z);
    }

    G4double getFermiMomentumConstant(const G4int /*A*/, const G4int /*Z*/) {
      return constantFermiMomentum;
    }

    G4double getFermiMomentumConstantLight(const G4int A, const G4int Z) {
// assert(Z>0 && A>0 && Z<=A);
      if(Z<clusterTableZSize && A<clusterTableASize) {
        const G4double rms = momentumRMS[Z][A];
        return ((rms>0.) ? rms : momentumRMS[6][12]) * Math::sqrtFiveThirds;
      } else
        return getFermiMomentumConstant(A,Z);
    }

    G4double getFermiMomentumMassDependent(const G4int A, const G4int /*Z*/) {
// assert(A>0);
      static const G4double alphaParam = 259.416; // MeV/c
      static const G4double betaParam  = 152.824; // MeV/c
      static const G4double gammaParam = 9.5157E-2;
      return alphaParam - betaParam*std::exp(-gammaParam*((G4double)A));
    }

    G4double getRPCorrelationCoefficient(const ParticleType t) {
// assert(t==Proton || t==Neutron || t==Lambda);
      return rpCorrelationCoefficient[t];
    }

    G4double getNeutronSkin() { return neutronSkin; }

    G4double getNeutronHalo() { return neutronHalo; }

    G4ThreadLocal G4double minDeltaMass = 0.;
    G4ThreadLocal G4double minDeltaMass2 = 0.;
    G4ThreadLocal G4double minDeltaMassRndm = 0.;
    G4ThreadLocal NuclearMassFn getTableMass = NULL;
    G4ThreadLocal ParticleMassFn getTableParticleMass = NULL;
    G4ThreadLocal SeparationEnergyFn getSeparationEnergy = NULL;
    G4ThreadLocal FermiMomentumFn getFermiMomentum = NULL;

    ParticleType getPionType(const G4int isosp) {
// assert(isosp == -2 || isosp == 0 || isosp == 2);
        if (isosp == -2) {
            return PiMinus;
        }
        else if (isosp == 0) {
            return PiZero;
        }
        else {
            return PiPlus;
        }
    }

    ParticleType getNucleonType(const G4int isosp) {
// assert(isosp == -1 || isosp == 1);
        if (isosp == -1) {
            return Neutron;
        }
        else {
            return Proton;
        }
    }

    ParticleType getDeltaType(const G4int isosp) {
// assert(isosp == -3 || isosp == -1 || isosp == 1 || isosp == 3);
        if (isosp == -3) {
            return DeltaMinus;
        }
        else if (isosp == -1) {
            return DeltaZero;
        }
        else if (isosp == 1) {
            return DeltaPlus;
        }
        else {
            return DeltaPlusPlus;
        }
    }


    ParticleType getSigmaType(const G4int isosp) {
// assert(isosp == -2 || isosp == 0 || isosp == 2);
        if (isosp == -2) {
            return SigmaMinus;
        }
        else if (isosp == 0) {
            return SigmaZero;
        }
        else {
            return SigmaPlus;
        }
    }

    ParticleType getKaonType(const G4int isosp) {
// assert(isosp == -1 || isosp == 1);
        if (isosp == -1) {
            return KZero;
        }
        else {
            return KPlus;
        }
    }

    ParticleType getAntiKaonType(const G4int isosp) {
// assert(isosp == -1 || isosp == 1);
        if (isosp == -1) {
            return KMinus;
        }
        else {
            return KZeroBar;
        }
    }

    G4double getWidth(const ParticleType pt) {
// assert(pt == PiPlus || pt == PiMinus || pt == PiZero || pt == Eta || pt == Omega || pt == EtaPrime || pt == KShort || pt == KLong || pt== KPlus || pt == KMinus || pt == Lambda || pt == SigmaPlus || pt == SigmaZero || pt == SigmaMinus);
          if(pt == PiPlus) {
              return piPlusWidth;
          } else if(pt == PiMinus) {
              return piMinusWidth;
          } else if(pt == PiZero) {
              return piZeroWidth;
          } else if(pt == Eta) {
              return etaWidth;
          } else if(pt == Omega) {
              return omegaWidth;
          } else if(pt == EtaPrime) {
              return etaPrimeWidth;
          } else if(pt == SigmaPlus) {
              return SigmaPlusWidth;
          } else if(pt == SigmaZero) {
              return SigmaZeroWidth;
          } else if(pt == SigmaMinus) {
              return SigmaMinusWidth;
          } else if(pt == KPlus) {
              return KPlusWidth;
          } else if(pt == KMinus) {
              return KMinusWidth;
          } else if(pt == KShort) {
              return KShortWidth;
          } else if(pt == KLong) {
              return KLongWidth;
          } else {
              INCL_ERROR("getWidth : Unknown particle type." << '\n');
              return 0.0;
          }
      }
      
  } // namespace ParticleTable
} // namespace G4INCL

