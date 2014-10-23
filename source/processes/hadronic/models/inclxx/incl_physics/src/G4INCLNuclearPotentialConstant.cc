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

/** \file G4INCLNuclearPotentialConstant.cc
 * \brief Isospin- and energy-independent nuclear potential.
 *
 * Provides a constant nuclear potential (V0).
 *
 * \date 17 January 2011
 * \author Davide Mancusi
 */

#include "G4INCLNuclearPotentialConstant.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  namespace NuclearPotential {

    // Constructors
    NuclearPotentialConstant::NuclearPotentialConstant(const G4int A, const G4int Z, const G4bool aPionPotential)
      : INuclearPotential(A, Z, aPionPotential)
    {
      initialize();
    }

    // Destructor
    NuclearPotentialConstant::~NuclearPotentialConstant() {
    }

    void NuclearPotentialConstant::initialize() {
      const G4double mp = ParticleTable::getINCLMass(Proton);
      const G4double mn = ParticleTable::getINCLMass(Neutron);

      const G4double theFermiMomentum = ParticleTable::getFermiMomentum(theA,theZ);

      fermiMomentum[Proton] = theFermiMomentum;
      const G4double theProtonFermiEnergy = std::sqrt(theFermiMomentum*theFermiMomentum + mp*mp) - mp;
      fermiEnergy[Proton] = theProtonFermiEnergy;

      fermiMomentum[Neutron] = theFermiMomentum;
      const G4double theNeutronFermiEnergy = std::sqrt(theFermiMomentum*theFermiMomentum + mn*mn) - mn;
      fermiEnergy[Neutron] = theNeutronFermiEnergy;

      fermiEnergy[DeltaPlusPlus] = fermiEnergy.find(Proton)->second;
      fermiEnergy[DeltaPlus] = fermiEnergy.find(Proton)->second;
      fermiEnergy[DeltaZero] = fermiEnergy.find(Neutron)->second;
      fermiEnergy[DeltaMinus] = fermiEnergy.find(Neutron)->second;

      const G4double theAverageSeparationEnergy = 0.5*(ParticleTable::getSeparationEnergy(Proton,theA,theZ)+ParticleTable::getSeparationEnergy(Neutron,theA,theZ));
      separationEnergy[Proton] = theAverageSeparationEnergy;
      separationEnergy[Neutron] = theAverageSeparationEnergy;

      // Use separation energies from the ParticleTable
      vNucleon = 0.5*(theProtonFermiEnergy + theNeutronFermiEnergy) + theAverageSeparationEnergy;
      vDelta = vNucleon;
      separationEnergy[DeltaPlusPlus] = vDelta - fermiEnergy.find(DeltaPlusPlus)->second;
      separationEnergy[DeltaPlus] = vDelta - fermiEnergy.find(DeltaPlus)->second;
      separationEnergy[DeltaZero] = vDelta - fermiEnergy.find(DeltaZero)->second;
      separationEnergy[DeltaMinus] = vDelta - fermiEnergy.find(DeltaMinus)->second;

      separationEnergy[PiPlus] = 0.;
      separationEnergy[PiZero] = 0.;
      separationEnergy[PiMinus] = 0.;

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

    G4double NuclearPotentialConstant::computePotentialEnergy(const Particle *particle) const {

      switch( particle->getType() )
      {
        case Proton:
        case Neutron:
          return vNucleon;
          break;

        case PiPlus:
        case PiZero:
        case PiMinus:
          return computePionPotentialEnergy(particle);
          break;

        case DeltaPlusPlus:
        case DeltaPlus:
        case DeltaZero:
        case DeltaMinus:
          return vDelta;
          break;
        case UnknownParticle:
          INCL_ERROR("Trying to compute potential energy of an unknown particle.");
          return 0.0;
          break;
        default:
          INCL_ERROR("Trying to compute potential energy of a malformed particle.");
          return 0.0;
          break;
      }
    }

  }
}

