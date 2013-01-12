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

#ifndef G4INCLCrossSections_hh
#define G4INCLCrossSections_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLIChannel.hh"

namespace G4INCL {
  class CrossSections {
  public:
    static G4double elastic(Particle const * const p1, Particle const * const p2);
    static G4double total(Particle const * const p1, Particle const * const p2);

    static G4double pionNucleon(Particle const * const p1, Particle const * const p2);

    static G4double recombination(Particle const * const p1, Particle const * const p2);
    static G4double deltaProduction(Particle const * const p1, Particle const * const p2);
    /** \brief Calculate the slope of the NN DDXS.
     *
     * \param energyCM energy in the CM frame, in MeV
     * \param iso total isospin of the system
     *
     * \return the slope of the angular distribution
     */
    static G4double calculateNNDiffCrossSection(G4double energyCM, G4int iso);

    /** \brief Compute the "interaction distance".
     *
     * Defined on the basis of the average value of the N-N cross sections at
     * the given kinetic energy.
     *
     * \return the interaction distance
     */
    static G4double interactionDistanceNN(const G4double projectileKineticEnergy);

    /** \brief Compute the "interaction distance".
     *
     * Defined on the basis of the average value of the pi-N cross sections at
     * the given kinetic energy.
     *
     * \return the interaction distance
     */
    static G4double interactionDistancePiN(const G4double projectileKineticEnergy);

    /** \brief The interaction distance for nucleons at 1 GeV.
     *
     * Used to determine the universe radius at any energy.
     */
    static G4double interactionDistanceNN1GeV() {
      static G4double answer = CrossSections::interactionDistanceNN(1000.);
      return answer;
    }

    /** \brief The interaction distance for pions at 1 GeV.
     *
     * Used to determine the universe radius at any energy.
     */
    static G4double interactionDistancePiN1GeV() {
      static G4double answer = CrossSections::interactionDistancePiN(1000.);
      return answer;
    }

  private:
    static G4double elasticNNHighEnergy(const G4double momentum);
    static G4double elasticProtonNeutron(const G4double momentum);
    static G4double elasticProtonProtonOrNeutronNeutron(const G4double momentum);
    static G4double elasticNN(Particle const * const p1, Particle const * const p2);
    static G4double elasticNNLegacy(Particle const * const p1, Particle const * const p2);

    static G4double deltaProduction(const G4int isospin, const G4double pCM);

    static G4double spnPiPlusPHE(const G4double x);
    static G4double spnPiMinusPHE(const G4double x);

  protected:
    CrossSections() {};
    ~CrossSections() {};
    
  };
}

#endif
