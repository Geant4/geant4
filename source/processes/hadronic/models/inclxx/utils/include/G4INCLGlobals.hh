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

#ifndef G4INCLGlobals_hh
#define G4INCLGlobals_hh 1

#include <cmath>
#include <string>
#include <vector>
#include "G4INCLParticleType.hh"

namespace G4INCL {
  class Particle;

  namespace PhysicalConstants {
    /// \brief \f$\hbar c\f$ [MeV*fm]
    const G4double hc = 197.328;

    /// \brief \f$\hbar^2 c^2\f$ [MeV^2*fm^2]
    const G4double hcSquared = hc*hc;

    /// \brief Fermi momentum [MeV/c]
    const G4double Pf = 1.37*hc;
    //  const G4double Pf = 1.36828*hc;

    /** \brief Coulomb conversion factor [MeV*fm]
     *
     * \f[ e^2/(4 pi epsilon_0) \f]
     */
    const G4double eSquared = 1.439964;
  }

  namespace Math {
    const G4double pi = 3.14159265358979323846264338328;
    const G4double twoPi = 2.0 * pi;
    const G4double tenPi = 10.0 * pi;
    const G4double piOverTwo = 0.5 * pi;
    const G4double oneOverSqrtTwo = 1./std::sqrt((G4double)2.);
    const G4double oneOverSqrtThree = 1./std::sqrt((G4double)3.);
    const G4double oneThird = 1./3.;
    const G4double twoThirds = 2./3.;
    const G4double sqrtFiveThirds = std::sqrt(5./3.);
    const G4double sqrtThreeFifths = std::sqrt(3./5.);

    inline G4double toDegrees(G4double radians) {
      return radians * (180.0 / pi);
    }

    inline G4int heaviside(G4int n) {
      if(n < 0) return 0;
      else return 1;
    }

    inline G4double pow13(G4double x) {
      return std::pow(x, oneThird);
    }

    inline G4double powMinus13(G4double x) {
      return std::pow(x, -oneThird);
    }

    inline G4double pow23(G4double x) {
      return std::pow(x, twoThirds);
    }

    inline G4double aSinH(G4double x) {
      return std::log(x + std::sqrt(x*x+1.));
    }

    /**
     * A simple sign function that allows us to port fortran code to c++ more easily.
     */
    template <typename T> inline G4int sign(const T t) {
      return t > 0 ? 1: t < 0 ? -1 : 0;
    }

    /// brief Return the largest of the two arguments
    template <typename T> inline T max(const T t1, const T t2) {
      return t1 > t2 ? t1 : t2;
    }

    /// brief Return the smallest of the two arguments
    template <typename T> inline T min(const T t1, const T t2) {
      return t1 < t2 ? t1 : t2;
    }

    /** \brief Cumulative distribution function for Gaussian
     *
     * A public-domain approximation taken from Abramowitz and Stegun. Applies
     * to a Gaussian with mean=0 and sigma=1.
     *
     * \param x a Gaussian variable
     */
    G4double gaussianCDF(const G4double x);

    /** \brief Generic cumulative distribution function for Gaussian
     *
     * A public-domain approximation taken from Abramowitz and Stegun. Applies
     * to a generic Gaussian.
     *
     * \param x a Gaussian variable
     * \param x0 mean of the Gaussian
     * \param sigma standard deviation of the Gaussian
     */
    G4double gaussianCDF(const G4double x, const G4double x0, const G4double sigma);

    /** \brief Inverse cumulative distribution function for Gaussian
     *
     * A public-domain approximation taken from Abramowitz and Stegun. Applies
     * to a Gaussian with mean=0 and sigma=1.
     *
     * \param x a uniform variate
     * \return a Gaussian variate
     */
    G4double inverseGaussianCDF(const G4double x);

    /// \brief Calculates arcsin with some tolerance on illegal arguments
    G4double arcSin(const G4double x);

    /// \brief Calculates arccos with some tolerance on illegal arguments
    G4double arcCos(const G4double x);
  }

  namespace ParticleConfig {
    G4bool isPair(Particle const * const p1, Particle const * const p2, ParticleType t1, ParticleType t2);
  }

#ifndef INCLXX_IN_GEANT4_MODE
  namespace String {
    void wrap(std::string &str, const size_t lineLength=78, const std::string &separators=" \t");
    void replaceAll(std::string &str, const std::string &from, const std::string &to, const size_t maxPosition=std::string::npos);
    std::vector<std::string> tokenize(std::string const &str, const std::string &delimiters);
    G4bool isInteger(std::string const &str);
    std::string expandPath(std::string const &path);
  }
#endif
}
#endif
