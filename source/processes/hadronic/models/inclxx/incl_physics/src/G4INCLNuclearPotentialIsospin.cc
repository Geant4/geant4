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

/** \file G4INCLNuclearPotentialIsospin.cc
 * \brief Isospin-dependent nuclear potential.
 *
 * Provides an isospin-dependent nuclear potential.
 *
 * Created on: 28 February 2011
 *     Author: Davide Mancusi
 */

#include "G4INCLNuclearPotentialIsospin.hh"
#include "G4INCLNuclearPotentialConstant.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  namespace NuclearPotential {

    // Constructors
    NuclearPotentialIsospin::NuclearPotentialIsospin(NuclearDensity *density, G4bool pionPotential)
      : INuclearPotential(density, pionPotential)
    {
      initialize();
    }

    // Destructor
    NuclearPotentialIsospin::~NuclearPotentialIsospin() {}

    void NuclearPotentialIsospin::initialize() {
      const G4double ZOverA = ((G4double) theDensity->getZ()) / ((G4double) theDensity->getA());

      const G4double mp = ParticleTable::getMass(Proton);
      fermiMomentum[Proton] = Pf * Math::pow13(2.*ZOverA);
      fermiEnergy[Proton] = std::sqrt(fermiMomentum[Proton]*fermiMomentum[Proton] + mp*mp) - mp;
      vProton = fermiEnergy[Proton] + ParticleTable::getSeparationEnergy(Proton);

      const G4double mn = ParticleTable::getMass(Neutron);
      fermiMomentum[Neutron] = Pf * Math::pow13(2.*(1.-ZOverA));
      fermiEnergy[Neutron] = std::sqrt(fermiMomentum[Neutron]*fermiMomentum[Neutron] + mn*mn) - mn;
      vNeutron = fermiEnergy[Neutron] + ParticleTable::getSeparationEnergy(Neutron);

      vDeltaPlus = vProton;
      vDeltaZero = vNeutron;
      vDeltaPlusPlus = 2*vDeltaPlus - vDeltaZero;
      vDeltaMinus = 2*vDeltaZero - vDeltaPlus;

      const G4double Tfpp = vDeltaPlusPlus - vProton + fermiEnergy.find(Proton)->second;
      const G4double Tfp = fermiEnergy.find(Proton)->second;
      const G4double Tf0 = fermiEnergy.find(Neutron)->second;
      const G4double Tfm = vDeltaMinus - vNeutron + fermiEnergy.find(Neutron)->second;
      fermiEnergy[DeltaPlusPlus] = Tfpp;
      fermiEnergy[DeltaPlus] = Tfp;
      fermiEnergy[DeltaZero] = Tf0;
      fermiEnergy[DeltaMinus] = Tfm;
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
	ERROR("No potential computed for particle of type Cluster.");
	return 0.0;
	break;
      case UnknownParticle:
	ERROR("Trying to compute potential energy for an unknown particle.");
	return 0.0;
	break;
      }

      ERROR("There is no potential for this type of particle.");
      return 0.0;
    }

  }
}

