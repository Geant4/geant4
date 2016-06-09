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

/** \file G4INCLINuclearPotential.hh
 * \brief Abstract G4interface to the nuclear potential.
 *
 * NuclearPotential-like classes should provide access to the value of the
 * potential of a particle in a particular context. For example, an instance of
 * a NuclearPotential class should be associated to every nucleus.
 *
 * Created on: 17 January 2011
 *     Author: Davide Mancusi
 */

#ifndef G4INCLINUCLEARPOTENTIAL_HH
#define G4INCLINUCLEARPOTENTIAL_HH 1

#include "G4INCLParticle.hh"
#include "G4INCLNuclearDensity.hh"
#include <map>
// #include <cassert>

namespace G4INCL {

  namespace NuclearPotential {

    class INuclearPotential {
      public:
        INuclearPotential(NuclearDensity * const nuclearDensity, G4bool pionPot)
          : theDensity(nuclearDensity), pionPotential(pionPot)
        {
          if(pionPotential) {
            const G4double ZOverA = ((G4double) theDensity->getZ()) / ((G4double) theDensity->getA());
            // As in INCL4.6, use the r0*A^(1/3) formula to estimate vc
            const G4double r = 1.12*Math::pow13((G4double)theDensity->getA());

            const G4double xsi = 1. - 2.*ZOverA;
            const G4double vc = 1.8*theDensity->getZ()/r; // 1.8 = 1.44*1.25
            vPiPlus = vPionDefault + 71.*xsi - vc;
            vPiZero = vPionDefault;
            vPiMinus = vPionDefault - 71.*xsi + vc;
          } else {
            vPiPlus = 0.0;
            vPiZero = 0.0;
            vPiMinus = 0.0;
          }

        }

        virtual ~INuclearPotential() {}

        inline NuclearDensity *getDensity() const {
          return theDensity;
        }

        void setDensity(NuclearDensity * const nuclearDensity) {
          theDensity = nuclearDensity;
        }

        /// \brief Do we have a pion potential?
        G4bool hasPionPotential() { return pionPotential; }

        virtual G4double computePotentialEnergy(const Particle * const p) const = 0;

        /** \brief Return the Fermi energy for a particle.
         *
         * \param p poG4inter to a Particle
         * \return Fermi energy for that particle type
         **/
        inline G4double getFermiEnergy(const Particle * const p) const { return fermiEnergy.find(p->getType())->second; }

        /** \brief Return the Fermi energy for a particle type.
         *
         * \param t particle type
         * \return Fermi energy for that particle type
         **/
        inline G4double getFermiEnergy(const ParticleType t) const { return fermiEnergy.find(t)->second; }

        /** \brief Return the Fermi momentum for a particle.
         *
         * \param p poG4inter to a Particle
         * \return Fermi momentum for that particle type
         **/
        inline G4double getFermiMomentum(const Particle * const p) const {
          if(p->isDelta()) {
            const G4double Tf = getFermiEnergy(p), m = p->getMass();
            return std::sqrt(Tf*(Tf+2.*m));
          } else
            return fermiMomentum.find(p->getType())->second;
        }

        /** \brief Return the Fermi momentum for a particle type.
         *
         * \param t particle type
         * \return Fermi momentum for that particle type
         **/
        inline G4double getFermiMomentum(const ParticleType t) const {
          // assert(t!=DeltaPlusPlus && t!=DeltaPlus && t!=DeltaZero && t!=DeltaMinus);
          return fermiMomentum.find(t)->second;
        }

      protected:
        /// \brief Compute the potential energy for the given pion.
        G4double computePionPotentialEnergy(const Particle * const p) const {
          // assert(p->getType()==PiPlus || p->getType()==PiZero || p->getType()==PiMinus);
          if(pionPotential && !p->isOutOfWell()) {
            switch( p->getType() ) {
              case PiPlus:
                return vPiPlus;
                break;
              case PiZero:
                return vPiZero;
                break;
              case PiMinus:
                return vPiMinus;
                break;
	    default: // Pion potential is defined and non-zero only for pions
	      return 0.0;
	      break;
            }
          }
          else
            return 0.0;
        }

        NuclearDensity *theDensity;

        /* \brief map of Fermi energies per particle type */
        std::map<ParticleType,G4double> fermiEnergy;
        /* \brief map of Fermi momenta per particle type */
        std::map<ParticleType,G4double> fermiMomentum;

      private:
        G4bool pionPotential;
        G4double vPiPlus, vPiZero, vPiMinus;
        static const G4double vPionDefault;

    };

  }

}

#endif /* G4INCLINUCLEARPOTENTIAL_HH_ */
