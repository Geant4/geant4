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


/** \file G4INCLINuclearPotential.cc
 * \brief Abstract interface to the nuclear potential.
 *
 * NuclearPotential-like classes should provide access to the value of the
 * potential of a particle in a particular context. For example, an instance of
 * a NuclearPotential class should be associated to every nucleus.
 *
 * \date 31 March 2011
 * \author Davide Mancusi
 */

#include "G4INCLINuclearPotential.hh"
#include "G4INCLNuclearPotentialEnergyIsospinSmooth.hh"
#include "G4INCLNuclearPotentialEnergyIsospin.hh"
#include "G4INCLNuclearPotentialIsospin.hh"
#include "G4INCLNuclearPotentialConstant.hh"

namespace G4INCL {

  namespace NuclearPotential {

    const G4double INuclearPotential::vPionDefault = 30.6; // MeV
    const G4double INuclearPotential::vKPlusDefault = -25.; // MeV (repulsive)
    const G4double INuclearPotential::vKMinusDefault = 60.; // MeV

    namespace {

      G4ThreadLocal std::map<long,INuclearPotential const *> *nuclearPotentialCache = NULL;

    }

    INuclearPotential const *createPotential(const PotentialType type, const G4int theA, const G4int theZ, const G4bool pionPotential) {
      if(!nuclearPotentialCache)
        nuclearPotentialCache = new std::map<long,INuclearPotential const *>;

      const long nuclideID = (pionPotential ? 1 : -1) * (1000*theZ + theA + 1000000*type); // MCNP-style nuclide IDs
      const std::map<long,INuclearPotential const *>::const_iterator mapEntry = nuclearPotentialCache->find(nuclideID);
      if(mapEntry == nuclearPotentialCache->end()) {
        INuclearPotential const *thePotential = NULL;
        switch(type) {
          case IsospinEnergySmoothPotential:
            thePotential = new NuclearPotentialEnergyIsospinSmooth(theA, theZ, pionPotential);
            break;
          case IsospinEnergyPotential:
            thePotential = new NuclearPotentialEnergyIsospin(theA, theZ, pionPotential);
            break;
          case IsospinPotential:
            thePotential = new NuclearPotentialIsospin(theA, theZ, pionPotential);
            break;
          case ConstantPotential:
            thePotential = new NuclearPotentialConstant(theA, theZ, pionPotential);
            break;
          default:
            INCL_FATAL("Unrecognized potential type at Nucleus creation." << '\n');
            break;
        }
        (*nuclearPotentialCache)[nuclideID] = thePotential;
        return thePotential;
      } else {
        return mapEntry->second;
      }
    }

    void clearCache() {
      if(nuclearPotentialCache) {
        for(std::map<long,INuclearPotential const *>::const_iterator i = nuclearPotentialCache->begin(), e=nuclearPotentialCache->end(); i!=e; ++i)
          delete i->second;
        nuclearPotentialCache->clear();
        delete nuclearPotentialCache;
        nuclearPotentialCache = NULL;
      }
    }

  } // namespace NuclearPotential

} // namespace G4INCL

