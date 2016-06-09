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

/** \file G4INCLNuclearPotentialConstant.cc
 * \brief Isospin- and energy-independent nuclear potential.
 *
 * Provides a constant nuclear potential (V0).
 *
 * Created on: 17 January 2011
 *     Author: Davide Mancusi
 */

#include "G4INCLNuclearPotentialConstant.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  namespace NuclearPotential {

    // Constructors
    NuclearPotentialConstant::NuclearPotentialConstant(NuclearDensity *density, G4bool pionPotential)
      : INuclearPotential(density, pionPotential)
    {
      initialize();
    }

    NuclearPotentialConstant::NuclearPotentialConstant(NuclearDensity *density, G4bool pionPotential, G4double /* nucleon */, G4double /* delta */)
      : INuclearPotential(density, pionPotential)
    {
      initialize();
    }

    // Destructor
    NuclearPotentialConstant::~NuclearPotentialConstant() {
    }

    void NuclearPotentialConstant::initialize() {
      const G4double mp = ParticleTable::getMass(Proton);
      fermiMomentum[Proton] = Pf;
      fermiEnergy[Proton] = std::sqrt(PfSquared + mp*mp) - mp;

      const G4double mn = ParticleTable::getMass(Neutron);
      fermiMomentum[Neutron] = Pf;
      fermiEnergy[Neutron] = std::sqrt(PfSquared + mn*mn) - mn;

      fermiEnergy[DeltaPlusPlus] = fermiEnergy.find(Proton)->second;
      fermiEnergy[DeltaPlus] = fermiEnergy.find(Proton)->second;
      fermiEnergy[DeltaZero] = fermiEnergy.find(Neutron)->second;
      fermiEnergy[DeltaMinus] = fermiEnergy.find(Neutron)->second;

      vNucleon = 0.5*(fermiEnergy[Proton]+fermiEnergy[Neutron])
        + 0.5*(ParticleTable::getSeparationEnergy(Proton)+ParticleTable::getSeparationEnergy(Neutron));
      vDelta = vNucleon;
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
          ERROR("Trying to compute potential energy of an unknown particle.");
          return 0.0;
          break;
        default:
          ERROR("Trying to compute potential energy of a malformed particle.");
          return 0.0;
          break;
      }
    }

  }
}

