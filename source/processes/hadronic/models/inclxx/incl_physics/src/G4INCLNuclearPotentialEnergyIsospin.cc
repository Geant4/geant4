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

/** \file G4INCLNuclearPotentialEnergyIsospin.cc
 * \brief Isospin- and energy-dependent nuclear potential.
 *
 * Provides an isospin- and energy-dependent nuclear potential.
 *
 * \date 21 March 2011
 * \author Davide Mancusi
 */

#include "G4INCLNuclearPotentialEnergyIsospin.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  namespace NuclearPotential {

    const G4double NuclearPotentialEnergyIsospin::alpha= 0.223;

    // Constructors
    NuclearPotentialEnergyIsospin::NuclearPotentialEnergyIsospin(const G4int A, const G4int Z, const G4bool aPionPotential)
      : NuclearPotentialIsospin(A,Z,aPionPotential)
    {}

    // Destructor
    NuclearPotentialEnergyIsospin::~NuclearPotentialEnergyIsospin() {}

    G4double NuclearPotentialEnergyIsospin::computePotentialEnergy(const Particle *particle) const {

      const G4double v0 = NuclearPotentialIsospin::computePotentialEnergy(particle);

      if(particle->isNucleon()) {
        const G4double t = particle->getKineticEnergy();
        const G4double tf = getFermiEnergy(particle);
        // Constant potential for T<Tf
        if(t < tf)
          return v0;

        // Linear function for T>Tf
        const G4double v = v0 - alpha*(t-tf)/(1-alpha);
        return (v>0.0) ? v : 0.0; // return 0.0 if v is negative
      } else
        return v0;
    }

  }
}

