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

/** \file G4INCLINuclearPotential.hh
 * \brief Abstract interface to the nuclear potential.
 *
 * NuclearPotential-like classes should provide access to the value of the
 * potential of a particle in a particular context. For example, an instance of
 * a NuclearPotential class should be associated to every nucleus.
 *
 * \date 17 January 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLINUCLEARPOTENTIAL_HH
#define G4INCLINUCLEARPOTENTIAL_HH 1

#include "G4INCLParticle.hh"
#include "G4INCLRandom.hh"
#include "G4INCLDeuteronDensity.hh"
#include <map>
// #include <cassert>

namespace G4INCL {

  namespace NuclearPotential {
    class INuclearPotential {
      public:
        INuclearPotential(const G4int A, const G4int Z, const G4bool pionPot) :
          theA(A),
          theZ(Z),
          pionPotential(pionPot)
        {
          if(pionPotential) {
            const G4double ZOverA = ((G4double) theZ) / ((G4double) theA);
            // As in INCL4.6, use the r0*A^(1/3) formula to estimate vc
            const G4double r = 1.12*Math::pow13((G4double)theA);

            const G4double xsi = 1. - 2.*ZOverA;
            const G4double vc = 1.25*PhysicalConstants::eSquared*theZ/r;
            vPiPlus = vPionDefault + 71.*xsi - vc;
            vPiZero = vPionDefault;
            vPiMinus = vPionDefault - 71.*xsi + vc;
            vKPlus = vKPlusDefault;
            vKZero = vKPlusDefault + 10.; // Hypothesis to be check
            vKMinus = vKMinusDefault;
            vKZeroBar = vKMinusDefault - 10.; // Hypothesis to be check
          } else {
            vPiPlus = 0.0;
            vPiZero = 0.0;
            vPiMinus = 0.0;
            vKPlus = 0.0;
            vKZero = 0.0;
            vKMinus = 0.0;
            vKZeroBar = 0.0;
          }
        }

        virtual ~INuclearPotential() {}

        /// \brief Do we have a pion potential?
        G4bool hasPionPotential() const { return pionPotential; }

        virtual G4double computePotentialEnergy(const Particle * const p) const = 0;

        /** \brief Return the Fermi energy for a particle.
         *
         * \param p pointer to a Particle
         * \return Fermi energy for that particle type
         **/
        inline G4double getFermiEnergy(const Particle * const p) const {
          std::map<ParticleType, G4double>::const_iterator i = fermiEnergy.find(p->getType());
// assert(i!=fermiEnergy.end());
          return i->second;
        }

        /** \brief Return the Fermi energy for a particle type.
         *
         * \param t particle type
         * \return Fermi energy for that particle type
         **/
        inline G4double getFermiEnergy(const ParticleType t) const {
          std::map<ParticleType, G4double>::const_iterator i = fermiEnergy.find(t);
// assert(i!=fermiEnergy.end());
          return i->second;
        }

        /** \brief Return the separation energy for a particle.
         *
         * \param p pointer to a Particle
         * \return separation energy for that particle type
         **/
        inline G4double getSeparationEnergy(const Particle * const p) const {
          std::map<ParticleType, G4double>::const_iterator i = separationEnergy.find(p->getType());
// assert(i!=separationEnergy.end());
          return i->second;
        }

        /** \brief Return the separation energy for a particle type.
         *
         * \param t particle type
         * \return separation energy for that particle type
         **/
        inline G4double getSeparationEnergy(const ParticleType t) const {
          std::map<ParticleType, G4double>::const_iterator i = separationEnergy.find(t);
// assert(i!=separationEnergy.end());
          return i->second;
        }

        /** \brief Return the Fermi momentum for a particle.
         *
         * \param p pointer to a Particle
         * \return Fermi momentum for that particle type
         **/
        inline G4double getFermiMomentum(const Particle * const p) const {
          if(p->isDelta()) {
            const G4double Tf = getFermiEnergy(p), mass = p->getMass();
            return std::sqrt(Tf*(Tf+2.*mass));
          } else {
            std::map<ParticleType, G4double>::const_iterator i = fermiMomentum.find(p->getType());
// assert(i!=fermiMomentum.end());
            return i->second;
          }
        }

        /** \brief Return the Fermi momentum for a particle type.
         *
         * \param t particle type
         * \return Fermi momentum for that particle type
         **/
        inline G4double getFermiMomentum(const ParticleType t) const {
// assert(t!=DeltaPlusPlus && t!=DeltaPlus && t!=DeltaZero && t!=DeltaMinus);
          std::map<ParticleType, G4double>::const_iterator i = fermiMomentum.find(t);
          return i->second;
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

      protected:
        /// \brief Compute the potential energy for the given kaon.
        G4double computeKaonPotentialEnergy(const Particle * const p) const {
// assert(p->getType()==KPlus || p->getType()==KZero || p->getType()==KZeroBar || p->getType()==KMinus|| p->getType()==KShort|| p->getType()==KLong);
          if(pionPotential && !p->isOutOfWell()) { // if pionPotental false -> kaonPotential false
            switch( p->getType() ) {
              case KPlus:
                return vKPlus;
                break;
              case KZero:
                return vKZero;
                break;
              case KZeroBar:
                return vKZeroBar;
                break;
              case KShort:
              case KLong:
                 return 0.0; // Should never be in the nucleus
                 break;
               case KMinus:
                 return vKMinus;
                 break;
               default:
                 return 0.0;
                 break;
               }
            }
        else
          return 0.0;
        }

      protected:
        /// \brief Compute the potential energy for the given pion resonances (Eta, Omega and EtaPrime and Gamma also).
        G4double computePionResonancePotentialEnergy(const Particle * const p) const {
// assert(p->getType()==Eta || p->getType()==Omega || p->getType()==EtaPrime || p->getType()==Photon);
          if(pionPotential && !p->isOutOfWell()) {
            switch( p->getType() ) {
              case Eta:
//jcd             return vPiZero;
//jcd              return vPiZero*1.5;
                return 0.0; // (JCD: seems to give better results)
                break;
              case Omega:
                return 15.0; // S.Friedrich et al., Physics Letters B736(2014)26-32. (V. Metag in Hyperfine Interact (2015) 234:25-31 gives 29 MeV)
                break;
              case EtaPrime:
                return 37.0; // V. Metag in Hyperfine Interact (2015) 234:25-31
                break;
              case Photon:
                return 0.0;
                break;
              default: 
                return 0.0;
                break;
            }
          }
          else
            return 0.0;
        }

      protected:
        /// \brief The mass number of the nucleus
        const G4int theA;
        /// \brief The charge number of the nucleus
        const G4int theZ;
      private:
        const G4bool pionPotential;
        G4double vPiPlus, vPiZero, vPiMinus;
        static const G4double vPionDefault;
        G4double vKPlus, vKZero, vKZeroBar, vKMinus;
        static const G4double vKPlusDefault;
        static const G4double vKMinusDefault;
      protected:
        /* \brief map of Fermi energies per particle type */
        std::map<ParticleType,G4double> fermiEnergy;
        /* \brief map of Fermi momenta per particle type */
        std::map<ParticleType,G4double> fermiMomentum;
        /* \brief map of separation energies per particle type */
        std::map<ParticleType,G4double> separationEnergy;

    };



    /** \brief Create an INuclearPotential object
     *
     * This is the method that should be used to instantiate objects derived
     * from INuclearPotential. It uses a caching mechanism to minimise
     * thrashing and speed up the code.
     *
     * \param type the type of the potential to be created
     * \param theA mass number of the nucleus
     * \param theZ charge number of the nucleus
     * \param pionPotential whether pions should also feel the potential
     * \return a pointer to the nuclear potential
     */
    INuclearPotential const *createPotential(const PotentialType type, const G4int theA, const G4int theZ, const G4bool pionPotential);

    /// \brief Clear the INuclearPotential cache
    void clearCache();

  }

}

#endif /* G4INCLINUCLEARPOTENTIAL_HH_ */
