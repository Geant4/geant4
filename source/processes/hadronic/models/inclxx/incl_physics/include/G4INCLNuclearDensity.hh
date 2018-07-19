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
#include "G4INCLInterpolationTable.hh"

namespace G4INCL {

  class NuclearDensity {
  public:
    NuclearDensity(const G4int A, const G4int Z, InterpolationTable const * const rpCorrelationTableProton, InterpolationTable const * const rpCorrelationTableNeutron);
    ~NuclearDensity();

    /// \brief Copy constructor
    NuclearDensity(const NuclearDensity &rhs);

    /// \brief Assignment operator
    NuclearDensity &operator=(const NuclearDensity &rhs);

    /// \brief Helper method for the assignment operator
    void swap(NuclearDensity &rhs);

    /** \brief Get the maximum allowed radius for a given momentum.
     *  \param t type of the particle
     *  \param p absolute value of the particle momentum, divided by the
     *           relevant Fermi momentum.
     *  \return maximum allowed radius.
     */
    G4double getMaxRFromP(const ParticleType t, const G4double p) const;

    G4double getMinPFromR(const ParticleType t, const G4double r) const;

    G4double getMaximumRadius() const { return theMaximumRadius; };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    G4double getTransmissionRadius(Particle const * const p) const {
      const ParticleType t = p->getType();
// assert(t!=Neutron && t!=PiZero && t!=DeltaZero && t!=Eta && t!=Omega && t!=EtaPrime && t!=Photon && t!= Lambda && t!=SigmaZero && t!=KZero && t!=KZeroBar && t!=KShort && t!=KLong); // no neutral particles here
      if(t==Composite) {
        return transmissionRadius[t] +
          ParticleTable::getNuclearRadius(t, p->getA(), p->getZ());
      } else
        return transmissionRadius[t];
    };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    G4double getTransmissionRadius(ParticleType type) const {
// assert(type!=Composite);
      return transmissionRadius[type];
    };

    /// \brief Get the mass number.
    G4int getA() const { return theA; }

    /// \brief Get the charge number.
    G4int getZ() const { return theZ; }

    G4double getProtonNuclearRadius() const { return theProtonNuclearRadius; }
    void setProtonNuclearRadius(const G4double r) { theProtonNuclearRadius = r; }

  private:

    /** \brief Initialize the transmission radius. */
    void initializeTransmissionRadii();

    G4int theA, theZ;
    G4double theMaximumRadius;
    /// \brief Represents INCL4.5's R0 variable
    G4double theProtonNuclearRadius;

    /* \brief map of transmission radii per particle type */
    G4double transmissionRadius[UnknownParticle];

    InterpolationTable const *rFromP[UnknownParticle];
    InterpolationTable const *pFromR[UnknownParticle];
  };

}

#endif
