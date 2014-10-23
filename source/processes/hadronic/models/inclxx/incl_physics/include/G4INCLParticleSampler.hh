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

/** \file G4INCLParticleSampler.hh
 * \brief Class for sampling particles in a nucleus
 *
 * \date 18 July 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLPARTICLESAMPLER_HH_
#define G4INCLPARTICLESAMPLER_HH_

#include "G4INCLNuclearDensity.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLInterpolationTable.hh"

namespace G4INCL {

  class ParticleSampler {

    public:
      /** \brief Constructor.
       *
       * \param A the mass number
       * \param Z the charge number
       */
      ParticleSampler(const G4int A, const G4int Z);

      /// \brief Destructor
      ~ParticleSampler();

      /// \brief Getter for theDensity
      NuclearDensity const *getDensity() const { return theDensity; }

      /// \brief Getter for thePotential
      NuclearPotential::INuclearPotential const *getPotential() const { return thePotential; }

      /// \brief Getter for rpCorrelationCoefficient
      G4double getRPCorrelationCoefficient(const ParticleType t) const {
// assert(t==Proton || t==Neutron);
        return rpCorrelationCoefficient[t];
      }

      /// \brief Setter for theDensity
      void setDensity(NuclearDensity const * const d);

      /// \brief Setter for thePotential
      void setPotential(NuclearPotential::INuclearPotential const * const p);

      /// \brief Setter for rpCorrelationCoefficient
      void setRPCorrelationCoefficient(const ParticleType t, const G4double corrCoeff) {
// assert(t==Proton || t==Neutron);
        rpCorrelationCoefficient[t] = corrCoeff;
      }

      ParticleList sampleParticles(ThreeVector const &position);
      void sampleParticlesIntoList(ThreeVector const &position, ParticleList &theList);

    private:

      void updateSampleOneParticleMethods();

      typedef Particle *(ParticleSampler::*ParticleSamplerMethod)(const ParticleType t) const;

      /** \brief Sample a list of particles.
       *
       * This method is a pointer to the method that does the real work for protons.
       */
      ParticleSamplerMethod sampleOneProton;

      /** \brief Sample a list of particles.
       *
       * This method is a pointer to the method that does the real work for neutrons.
       */
      ParticleSamplerMethod sampleOneNeutron;

      /// \brief Sample one particle taking into account the rp-correlation
      Particle *sampleOneParticleWithRPCorrelation(const ParticleType t) const;

      /// \brief Sample one particle not taking into account the rp-correlation
      Particle *sampleOneParticleWithoutRPCorrelation(const ParticleType t) const;

      /// \brief Sample one particle with a fuzzy rp-correlation
      Particle *sampleOneParticleWithFuzzyRPCorrelation(const ParticleType t) const;

      /// \brief Mass number
      const G4int theA;

      /// \brief Charge number
      const G4int theZ;

      /// \brief Array of pointers to the r-space CDF table
      InterpolationTable const *theRCDFTable[UnknownParticle];

      /// \brief Array of pointers to the p-space CDF table
      InterpolationTable const *thePCDFTable[UnknownParticle];

      /// \brief Pointer to the Cluster's NuclearDensity
      NuclearDensity const *theDensity;

      /// \brief Pointer to the Cluster's NuclearPotential
      NuclearPotential::INuclearPotential const *thePotential;

      /// \brief Correlation coefficients for the r-p correlation
      G4double rpCorrelationCoefficient[UnknownParticle];
  };

}

#endif // G4INCLPARTICLESAMPLER_HH_
