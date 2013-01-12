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

namespace G4INCL {

  class ParticleSampler {

    public:
      /** \brief Constructor.
       *
       * \param A the mass number
       * \param Z the charge number
       * \param rCDFTable the interpolation table for the inverse CDF in r-space
       * \param pCDFTable the interpolation table for the inverse CDF in p-space
       */
      ParticleSampler(const G4int A, const G4int Z, InverseInterpolationTable const * const rCDFTable, InverseInterpolationTable const * const pCDFTable);

      /// \brief Destructor
      ~ParticleSampler();

      /// \brief Getter for theDensity
      NuclearDensity const *getDensity() const { return theDensity; }

      /// \brief Getter for thePotential
      NuclearPotential::INuclearPotential const *getPotential() const { return thePotential; }

      /// \brief Setter for theDensity
      void setDensity(NuclearDensity const * const d);

      /// \brief Setter for thePotential
      void setPotential(NuclearPotential::INuclearPotential const * const p);

      ParticleList sampleParticles(ThreeVector const &position) const;

    private:

      void updateSampleOneParticleMethod();

      typedef Particle *(ParticleSampler::*ParticleSamplerMethod)(const ParticleType t) const;
      /** \brief Sample a list of particles.
       *
       * This method is a pointer to the method that does the real work.
       */
      ParticleSamplerMethod sampleOneParticle;

      /// \brief Sample one particle taking into account the rp-correlation
      Particle *sampleOneParticleWithRPCorrelation(const ParticleType t) const;

      /// \brief Sample one particle not taking into account the rp-correlation
      Particle *sampleOneParticleWithoutRPCorrelation(const ParticleType t) const;

      /// \brief Mass number
      const G4int theA;

      /// \brief Charge number
      const G4int theZ;

      /// \brief Pointer to the r-space CDF table
      InverseInterpolationTable const *theRCDFTable;

      /// \brief Pointer to the p-space CDF table
      InverseInterpolationTable const *thePCDFTable;

      /// \brief Pointer to the Cluster's NuclearDensity
      NuclearDensity const *theDensity;

      /// \brief Pointer to the Cluster's NuclearPotential
      NuclearPotential::INuclearPotential const *thePotential;
  };

}

#endif // G4INCLPARTICLESAMPLER_HH_
