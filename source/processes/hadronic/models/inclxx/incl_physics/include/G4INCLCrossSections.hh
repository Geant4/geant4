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

#ifndef G4INCLCrossSections_hh
#define G4INCLCrossSections_hh 1

#include "G4INCLICrossSections.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {
  namespace CrossSections {
      G4double elastic(Particle const * const p1, Particle const * const p2);
      G4double total(Particle const * const p1, Particle const * const p2);

      G4double NDeltaToNN(Particle const * const p1, Particle const * const p2);
      G4double NNToNDelta(Particle const * const p1, Particle const * const p2);
      G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2);
      G4double piNToDelta(Particle const * const p1, Particle const * const p2);
      G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2);
      G4double piNToEtaN(Particle const * const p1, Particle const * const p2);
      G4double piNToOmegaN(Particle const * const p1, Particle const * const p2);
      G4double piNToEtaPrimeN(Particle const * const p1, Particle const * const p2);
		   	G4double etaNToPiN(Particle const * const p1, Particle const * const p2);
		   	G4double etaNToPiPiN(Particle const * const p1, Particle const * const p2);
      G4double omegaNToPiN(Particle const * const p1, Particle const * const p2);
      G4double omegaNToPiPiN(Particle const * const p1, Particle const * const p2);
      G4double etaPrimeNToPiN(Particle const * const p1, Particle const * const p2);

      G4double NNToNNEta(Particle const * const p1, Particle const * const p2);
      G4double NNToNNEtaExclu(Particle const * const p1, Particle const * const p2);
      G4double NNToNNEtaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2);
      G4double NNToNDeltaEta(Particle const * const p1, Particle const * const p2);
      G4double NNToNNOmega(Particle const * const p1, Particle const * const p2);
      G4double NNToNNOmegaExclu(Particle const * const p1, Particle const * const p2);
      G4double NNToNNOmegaxPi(const G4int xpi, Particle const * const p1, Particle const * const p2);
      G4double NNToNDeltaOmega(Particle const * const p1, Particle const * const p2);
      
      /// \brief Strange cross sections
      G4double NNToNLK(Particle const * const p1, Particle const * const p2);
      G4double NNToNSK(Particle const * const p1, Particle const * const p2);
      G4double NNToNLKpi(Particle const * const p1, Particle const * const p2);
      G4double NNToNSKpi(Particle const * const p1, Particle const * const p2);
      G4double NNToNLK2pi(Particle const * const p1, Particle const * const p2);
      G4double NNToNSK2pi(Particle const * const p1, Particle const * const p2);
      G4double NNToNNKKb(Particle const * const p1, Particle const * const p2);
      G4double NNToMissingStrangeness(Particle const * const p1, Particle const * const p2);
		G4double NDeltaToNLK(Particle const * const p1, Particle const * const p2);
		G4double NDeltaToNSK(Particle const * const p1, Particle const * const p2);
		G4double NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2);
		G4double NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2);
		G4double NDeltaToNNKKb(Particle const * const p1, Particle const * const p2);
      G4double NpiToLK(Particle const * const p1, Particle const * const p2);
      G4double NpiToSK(Particle const * const p1, Particle const * const p2);
		G4double p_pimToSzKz(Particle const * const p1, Particle const * const p2);
		G4double p_pimToSmKp(Particle const * const p1, Particle const * const p2);
		G4double p_pizToSzKp(Particle const * const p1, Particle const * const p2);
      G4double NpiToLKpi(Particle const * const p1, Particle const * const p2);
      G4double NpiToSKpi(Particle const * const p1, Particle const * const p2);
      G4double NpiToLK2pi(Particle const * const p1, Particle const * const p2);
      G4double NpiToSK2pi(Particle const * const p1, Particle const * const p2);
      G4double NpiToNKKb(Particle const * const p1, Particle const * const p2);
      G4double NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2);
         G4double NLToNS(Particle const * const p1, Particle const * const p2);
         G4double NSToNL(Particle const * const p1, Particle const * const p2);
         G4double NSToNS(Particle const * const p1, Particle const * const p2);
      G4double NKToNK(Particle const * const p1, Particle const * const p2);
      G4double NKToNKpi(Particle const * const p1, Particle const * const p2);
      G4double NKToNK2pi(Particle const * const p1, Particle const * const p2);
         G4double NKbToNKb(Particle const * const p1, Particle const * const p2);
         G4double NKbToSpi(Particle const * const p1, Particle const * const p2);
         G4double NKbToLpi(Particle const * const p1, Particle const * const p2);
         G4double NKbToS2pi(Particle const * const p1, Particle const * const p2);
         G4double NKbToL2pi(Particle const * const p1, Particle const * const p2);
         G4double NKbToNKbpi(Particle const * const p1, Particle const * const p2);
         G4double NKbToNKb2pi(Particle const * const p1, Particle const * const p2);
      G4double NYelastic(Particle const * const p1, Particle const * const p2);
      G4double NKbelastic(Particle const * const p1, Particle const * const p2);
      G4double NKelastic(Particle const * const p1, Particle const * const p2);
      
      /** \brief Calculate the slope of the NN DDXS.
       *
       * \param energyCM energy in the CM frame, in MeV
       * \param iso total isospin of the system
       *
       * \return the slope of the angular distribution, in (GeV/c)^(-2)
       */
      G4double calculateNNAngularSlope(G4double energyCM, G4int iso);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the N-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      G4double interactionDistanceNN(const ParticleSpecies &aSpecies, const G4double kineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the pi-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      G4double interactionDistancePiN(const G4double projectileKineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the K-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      G4double interactionDistanceKN(const G4double projectileKineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the Kbar-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      G4double interactionDistanceKbarN(const G4double projectileKineticEnergy);

      /** \brief Compute the "interaction distance".
       *
       * Defined on the basis of the average value of the Y-N cross sections at
       * the given kinetic energy.
       *
       * \return the interaction distance
       */
      G4double interactionDistanceYN(const G4double projectileKineticEnergy);

      void setCrossSections(ICrossSections *c);

      void deleteCrossSections();

      void initialize(Config const * const theConfig);

  }
}

#endif
