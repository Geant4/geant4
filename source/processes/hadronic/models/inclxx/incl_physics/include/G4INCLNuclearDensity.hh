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

#ifndef G4INCLNuclearDensity_hh
#define G4INCLNuclearDensity_hh 1

#include <vector>
#include <map>
// #include <cassert>
#include "G4INCLThreeVector.hh"
#include "G4INCLIFunction1D.hh"
#include "G4INCLParticle.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLRandom.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLInverseInterpolationTable.hh"

namespace G4INCL {

  class NuclearDensity {
  public:
    NuclearDensity(G4int A, G4int Z, InverseInterpolationTable *rpCorrelationTable);
    ~NuclearDensity();

    /// \brief Copy constructor
    NuclearDensity(const NuclearDensity &rhs);

    /// \brief Assignment operator
    NuclearDensity &operator=(const NuclearDensity &rhs);

    /// \brief Helper method for the assignment operator
    void swap(NuclearDensity &rhs);

    /** \brief Get the maximum allowed radius for a given momentum.
     *  \param p Absolute value of the particle momentum, divided by the
     *  relevant Fermi momentum.
     *  \return Maximum allowed radius.
     */
    G4double getMaxRFromP(G4double p) const;

    G4double getMaxTFromR(G4double r) const;

    G4double getMaximumRadius() const { return theMaximumRadius; };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    G4double getTransmissionRadius(Particle const * const p) const {
      const ParticleType t = p->getType();
// assert(t!=Neutron && t!=PiZero && t!=DeltaZero); // no neutral particles here
      if(t==Composite) {
        return transmissionRadius[t] +
          ParticleTable::getNuclearRadius(p->getA(), p->getZ());
      } else
        return transmissionRadius[t];
    };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    G4double getTransmissionRadius(ParticleType type) {
// assert(type!=Composite);
      return transmissionRadius[type];
    };

    /// \brief Get the mass number.
    G4int getA() const { return theA; }

    /// \brief Get the charge number.
    G4int getZ() const { return theZ; }

    G4double getNuclearRadius() { return theNuclearRadius; }

  private:

    /** \brief Initialize the transmission radius. */
    void initializeTransmissionRadii();

    G4int theA, theZ;
    G4double theMaximumRadius;
    /// \brief Represents INCL4.5's R0 variable
    G4double theNuclearRadius;

    /* \brief map of transmission radii per particle type */
    G4double transmissionRadius[UnknownParticle];

    InverseInterpolationTable *rFromP;
    InverseInterpolationTable *tFromR;
  };

}

#endif
