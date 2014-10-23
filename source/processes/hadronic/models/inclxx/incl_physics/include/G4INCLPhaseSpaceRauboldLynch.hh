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

#ifndef G4INCLPHASESPACERAUBOLDLYNCH_HH
#define G4INCLPHASESPACERAUBOLDLYNCH_HH

#include "G4INCLThreeVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIPhaseSpaceGenerator.hh"
#include "G4INCLInterpolationTable.hh"
#include <vector>

namespace G4INCL {
  /// \brief Generate momenta using the RauboldLynch method
  class PhaseSpaceRauboldLynch : public IPhaseSpaceGenerator {

    public:
      PhaseSpaceRauboldLynch();
      virtual ~PhaseSpaceRauboldLynch();

      /// \brief Dummy copy constructor to silence Coverity warning
      PhaseSpaceRauboldLynch(PhaseSpaceRauboldLynch const &other);

      /// \brief Dummy assignment operator to silence Coverity warning
      PhaseSpaceRauboldLynch &operator=(PhaseSpaceRauboldLynch const &rhs);

      /** \brief Generate momenta according to a uniform, Lorentz-invariant phase-space model
       *
       * This function will assign momenta to the particles in the list that is
       * passed as an argument. The event is generated in the CM frame.
       *
       * \param sqrtS total centre-of-mass energy of the system
       * \param particles list of particles
       */
      void generate(const G4double sqrtS, ParticleList &particles);

      /// \brief Return the largest generated weight
      G4double getMaxGeneratedWeight() const;

    private:
      std::vector<G4double> masses;
      std::vector<G4double> sumMasses;
      std::vector<G4double> rnd;
      std::vector<G4double> invariantMasses;
      std::vector<G4double> momentaCM;
      size_t nParticles;
      G4double sqrtS;
      G4double availableEnergy;
      G4double maxGeneratedWeight;

      static const size_t nMasslessParticlesTable = 13;
      static const size_t wMaxNE = 30;
      static const G4double wMaxMasslessX[wMaxNE];
      static const G4double wMaxMasslessY[wMaxNE];
      static const G4double wMaxCorrectionX[wMaxNE];
      static const G4double wMaxCorrectionY[wMaxNE];
      static const G4double wMaxInterpolationMargin;

      InterpolationTable *wMaxMassless;
      InterpolationTable *wMaxCorrection;

      static const size_t wMaxNP = 20;

      /// \brief Precalculated coefficients: -ln(n)
      G4double prelog[wMaxNP];

      /// \brief Initialize internal structures (masses and sum of masses)
      void initialize(ParticleList &particles);

      /// \brief Compute the maximum possible weight using a naive algorithm
      G4double computeMaximumWeightNaive();

      /// \brief Compute the maximum possible weight using parametrizations
      G4double computeMaximumWeightParam();

      /// \brief Compute the maximum possible weight
      G4double computeWeight();

      /// \brief Generate an event
      void generateEvent(ParticleList &particles);
  };
}

#endif // G4INCLPHASESPACERAUBOLDLYNCH_HH
