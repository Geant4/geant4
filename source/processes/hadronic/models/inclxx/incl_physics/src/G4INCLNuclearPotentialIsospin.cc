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

/** \file G4INCLNuclearPotentialIsospin.cc
 * \brief Isospin-dependent nuclear potential.
 *
 * Provides an isospin-dependent nuclear potential.
 *
 * \date 28 February 2011
 * \author Davide Mancusi
 */

#include "G4INCLNuclearPotentialIsospin.hh"
#include "G4INCLNuclearPotentialConstant.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  namespace NuclearPotential {

    // Constructors
    NuclearPotentialIsospin::NuclearPotentialIsospin(const G4int A, const G4int Z, const G4bool aPionPotential)
      : INuclearPotential(A, Z, aPionPotential)
    {
      initialize();
    }

    // Destructor
    NuclearPotentialIsospin::~NuclearPotentialIsospin() {}

    void NuclearPotentialIsospin::initialize() {
      const G4double ZOverA = ((G4double) theZ) / ((G4double) theA);

      const G4double mp = ParticleTable::getINCLMass(Proton);
      const G4double mn = ParticleTable::getINCLMass(Neutron);

      const G4double theFermiMomentum = ParticleTable::getFermiMomentum(theA,theZ);

      fermiMomentum[Proton] = theFermiMomentum * Math::pow13(2.*ZOverA);
      const G4double theProtonFermiEnergy = std::sqrt(fermiMomentum[Proton]*fermiMomentum[Proton] + mp*mp) - mp;
      fermiEnergy[Proton] = theProtonFermiEnergy;
      // Use separation energies from the ParticleTable
      const G4double theProtonSeparationEnergy = ParticleTable::getSeparationEnergy(Proton,theA,theZ);
      separationEnergy[Proton] = theProtonSeparationEnergy;
      vProton = theProtonFermiEnergy + theProtonSeparationEnergy;

      fermiMomentum[Neutron] = theFermiMomentum * Math::pow13(2.*(1.-ZOverA));
      const G4double theNeutronFermiEnergy = std::sqrt(fermiMomentum[Neutron]*fermiMomentum[Neutron] + mn*mn) - mn;
      fermiEnergy[Neutron] = theNeutronFermiEnergy;
      // Use separation energies from the ParticleTable
      const G4double theNeutronSeparationEnergy = ParticleTable::getSeparationEnergy(Neutron,theA,theZ);
      separationEnergy[Neutron] = theNeutronSeparationEnergy;
      vNeutron = theNeutronFermiEnergy + theNeutronSeparationEnergy;

      const G4double separationEnergyDeltaPlusPlus = 2.*theProtonSeparationEnergy - theNeutronSeparationEnergy;
      separationEnergy[DeltaPlusPlus] = separationEnergyDeltaPlusPlus;
      separationEnergy[DeltaPlus] = theProtonSeparationEnergy;
      separationEnergy[DeltaZero] = theNeutronSeparationEnergy;
      const G4double separationEnergyDeltaMinus = 2.*theNeutronSeparationEnergy - theProtonSeparationEnergy;
      separationEnergy[DeltaMinus] = separationEnergyDeltaMinus;

      const G4double tinyMargin = 1E-7;
      vDeltaPlus = vProton;
      vDeltaZero = vNeutron;
      vDeltaPlusPlus = std::max(separationEnergyDeltaPlusPlus + tinyMargin, 2.*vDeltaPlus - vDeltaZero);
      vDeltaMinus = std::max(separationEnergyDeltaMinus + tinyMargin, 2.*vDeltaZero - vDeltaPlus);

      separationEnergy[PiPlus] = theProtonSeparationEnergy - theNeutronSeparationEnergy;
      separationEnergy[PiZero] = 0.;
      separationEnergy[PiMinus] = theNeutronSeparationEnergy - theProtonSeparationEnergy;

      separationEnergy[Eta]      = 0.;
      separationEnergy[Omega]    = 0.;
      separationEnergy[EtaPrime] = 0.;
      separationEnergy[Photon]   = 0.;

      fermiEnergy[DeltaPlusPlus] = vDeltaPlusPlus - separationEnergy[DeltaPlusPlus];
      fermiEnergy[DeltaPlus] = vDeltaPlus - separationEnergy[DeltaPlus];
      fermiEnergy[DeltaZero] = vDeltaZero - separationEnergy[DeltaZero];
      fermiEnergy[DeltaMinus] = vDeltaMinus - separationEnergy[DeltaMinus];

      INCL_DEBUG("Table of separation energies [MeV] for A=" << theA << ", Z=" << theZ << ":" << '\n'
            << "  proton:  " << separationEnergy[Proton] << '\n'
            << "  neutron: " << separationEnergy[Neutron] << '\n'
            << "  delta++: " << separationEnergy[DeltaPlusPlus] << '\n'
            << "  delta+:  " << separationEnergy[DeltaPlus] << '\n'
            << "  delta0:  " << separationEnergy[DeltaZero] << '\n'
            << "  delta-:  " << separationEnergy[DeltaMinus] << '\n'
            << "  pi+:     " << separationEnergy[PiPlus] << '\n'
            << "  pi0:     " << separationEnergy[PiZero] << '\n'
            << "  pi-:     " << separationEnergy[PiMinus] << '\n'
            << "  eta:     " << separationEnergy[Eta] << '\n'
            << "  omega:   " << separationEnergy[Omega] << '\n'
            << "  etaprime:" << separationEnergy[EtaPrime] << '\n'
            << "  photon:  " << separationEnergy[Photon] << '\n'
            );

      INCL_DEBUG("Table of Fermi energies [MeV] for A=" << theA << ", Z=" << theZ << ":" << '\n'
            << "  proton:  " << fermiEnergy[Proton] << '\n'
            << "  neutron: " << fermiEnergy[Neutron] << '\n'
            << "  delta++: " << fermiEnergy[DeltaPlusPlus] << '\n'
            << "  delta+:  " << fermiEnergy[DeltaPlus] << '\n'
            << "  delta0:  " << fermiEnergy[DeltaZero] << '\n'
            << "  delta-:  " << fermiEnergy[DeltaMinus] << '\n'
            );

      INCL_DEBUG("Table of Fermi momenta [MeV/c] for A=" << theA << ", Z=" << theZ << ":" << '\n'
            << "  proton:  " << fermiMomentum[Proton] << '\n'
            << "  neutron: " << fermiMomentum[Neutron] << '\n'
            );
    }

    G4double NuclearPotentialIsospin::computePotentialEnergy(const Particle *particle) const {

      switch( particle->getType() )
      {
        case Proton:
          return vProton;
          break;
        case Neutron:
          return vNeutron;
          break;

        case PiPlus:
        case PiZero:
        case PiMinus:
          return computePionPotentialEnergy(particle);
          break;

        case Eta:
        case Omega:
		case EtaPrime:
          return computePionResonancePotentialEnergy(particle);
          break;

        case Photon:
          return 0.0;
          break;

        case DeltaPlusPlus:
          return vDeltaPlusPlus;
          break;
        case DeltaPlus:
          return vDeltaPlus;
          break;
        case DeltaZero:
          return vDeltaZero;
          break;
        case DeltaMinus:
          return vDeltaMinus;
          break;
      case Composite:
	INCL_ERROR("No potential computed for particle of type Cluster.");
	return 0.0;
	break;
      case UnknownParticle:
	INCL_ERROR("Trying to compute potential energy for an unknown particle.");
	return 0.0;
	break;
      }

      INCL_ERROR("There is no potential for this type of particle.");
      return 0.0;
    }

  }
}

