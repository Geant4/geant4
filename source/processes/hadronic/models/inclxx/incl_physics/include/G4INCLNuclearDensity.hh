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

#ifndef G4INCLNuclearDensity_hh
#define G4INCLNuclearDensity_hh 1

#include <vector>
#include <map>
#include "G4INCLThreeVector.hh"
#include "G4INCLIFunction.hh"
#include "G4INCLParticle.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  class NuclearDensity {
  public:
    NuclearDensity(G4int A, G4int Z, IFunction1D *densityFunction);
    NuclearDensity(G4int A, G4int Z, IFunction1D *densityFunction,
		     G4double radius, G4double maxRadius, G4double diffuseness);
    // NuclearDensity(G4int A, G4int Z, IFunction1D *densityFunction,
    // 		     G4double radius, G4double maxRadius, G4double diffuseness);
    ~NuclearDensity();

    G4double getFirstDerivative(G4int index) const;

    /** \brief Get the maximum allowed radius for a given momentum.
     *  \param p Absolute value of the particle momentum, divided by the
     *  relevant Fermi momentum.
     *  \return Maximum allowed radius.
     */
    G4double getMaxRFromP(G4double p) const;
    G4double getMaxRFromPLegacy(G4double p) const;
    G4double getMaxRFromPNew(G4double p) const;

    G4double getMaxTFromR(G4double r) const;

    G4double getMaximumRadius() const { return theMaximumRadius; };

    /** \brief Initialize the transmission radius. */
    void initializeTransmissionRadii();

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    G4double getTransmissionRadius(Particle const * const p) const {
      if(p->getType()==Composite) {
        return transmissionRadius.find(p->getType())->second +
          ParticleTable::getClusterRMS(p->getA(), p->getZ());
      } else
        return transmissionRadius.find(p->getType())->second;
    };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    G4double getTransmissionRadius(ParticleType type) {
      return transmissionRadius[type];
    };

    /// \brief Get the mass number.
    G4int getA() const { return theA; }

    /// \brief Get the charge number.
    G4int getZ() const { return theZ; }

    G4double getCentralRadius() { return theCentralRadius; }

  private:
    /**
     * New implementation of the density G4interpolation function
     * without gotos.
     */
    G4double getDensityNew(G4double) const;

    /**
     * Direct translation of the FORTRAN version of the density
     * G4interpolation routine.
     */
    G4double getDensityLegacy(G4double) const;

    void initializeDensity();
    void initializeFirstDerivative();
    G4double G4integrate(G4double ami, G4double ama, G4double step) const;
    void initMaterial(G4int iamat, G4int izmat);

    G4int theA, theZ;
    IFunction1D *densityFunction;
    G4double theRadiusParameter, theMaximumRadius, theDiffusenessParameter;
    /// \brief Represents INCL4.5's R0 variable
    G4double theCentralRadius;

    void computeCentralRadius() {
      if(theA>=6 && theA<19)
        theCentralRadius = 1.581*theDiffusenessParameter*
          (2.+5.*theRadiusParameter)/(2.+3.*theRadiusParameter);
      else
        theCentralRadius = theRadiusParameter;
    }

    /* \brief map of transmission radii per particle type */
    std::map<ParticleType,G4double> transmissionRadius;

    std::vector<G4double> x, y, s;
    std::vector<G4double> r_t, tmin, s_loce;
  };

}

#endif
