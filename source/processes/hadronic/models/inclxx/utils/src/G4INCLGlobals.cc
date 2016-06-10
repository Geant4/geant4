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

#include "G4INCLGlobals.hh"
#include "G4INCLParticle.hh"
#ifdef HAVE_WORDEXP
#include <wordexp.h>
#endif

namespace G4INCL {
  namespace Math {

    namespace {

      // constants for the Gaussian CDF approximation
      const G4double gcdfa1 =  0.254829592;
      const G4double gcdfa2 = -0.284496736;
      const G4double gcdfa3 =  1.421413741;
      const G4double gcdfa4 = -1.453152027;
      const G4double gcdfa5 =  1.061405429;
      const G4double gcdfp  =  0.3275911;

      // constants for the inverse Gaussian CDF approximation
      const G4double igcdfc1 = 2.515517;
      const G4double igcdfc2 = 0.802853;
      const G4double igcdfc3 = 0.010328;
      const G4double igcdfd1 = 1.432788;
      const G4double igcdfd2 = 0.189269;
      const G4double igcdfd3 = 0.001308;

      G4double inverseGaussianCDFRational(const G4double t) {
        // Abramowitz and Stegun formula 26.2.23.
        // The absolute value of the error should be less than 4.5 e-4.
        return t - ((igcdfc3*t + igcdfc2)*t + igcdfc1) /
          (((igcdfd3*t + igcdfd2)*t + igcdfd1)*t + 1.0);
      }

    }

    G4double gaussianCDF(const G4double x)
    {
      // Save the sign of x
      const G4double sgn = sign(x);
      const G4double z = std::fabs(x) * oneOverSqrtTwo;

      // A&S formula 7.1.26
      G4double t = 1.0/(1.0 + gcdfp*z);
      G4double y = 1.0 - (((((gcdfa5*t + gcdfa4)*t) + gcdfa3)*t + gcdfa2)*t + gcdfa1)*t*std::exp(-z*z);

      return 0.5*(1.0 + sgn*y);
    }

    G4double gaussianCDF(const G4double x, const G4double x0, const G4double sigma) {
      return gaussianCDF((x-x0)/sigma);
    }

    G4double inverseGaussianCDF(G4double x) {
      if (x < 0.5)
        return -inverseGaussianCDFRational( std::sqrt(-2.0*std::log(x)) );
      else
        return inverseGaussianCDFRational( std::sqrt(-2.0*std::log(1.-x)) );
    }

    G4double arcSin(const G4double x) {
// assert(x>-1.000001 && x<1.000001);
      return ((x > 1.) ? 0. : ((x<-1.) ? pi : std::asin(x)));
    }

    G4double arcCos(const G4double x) {
// assert(x>-1.000001 && x<1.000001);
      return ((x > 1.) ? 0. : ((x<-1.) ? pi : std::acos(x)));
    }
  }

  namespace ParticleConfig {
    G4bool isPair(Particle const * const p1, Particle const * const p2, ParticleType t1, ParticleType t2) {
      return ((p1->getType() == t1 && p2->getType() == t2) ||
	      (p1->getType() == t2 && p2->getType() == t1));
    }
  }

#ifndef INCLXX_IN_GEANT4_MODE
  namespace String {
    void wrap(std::string &str, const size_t lineLength, const std::string &separators) {
      const size_t len = str.size();
      size_t startPos = 0;
      while(len-startPos > lineLength) { /* Loop checking, 10.07.2015, D.Mancusi */
        const size_t nextNewline = str.find('\n', startPos);
        if(nextNewline!=std::string::npos && nextNewline-startPos<=lineLength)
          startPos = nextNewline+1;
        else {
          size_t lastSeparator = str.find_last_of(separators, startPos+lineLength);
          if(lastSeparator!=std::string::npos)
            str[lastSeparator] = '\n';
          startPos = lastSeparator+1;
        }
      }
    }

    void replaceAll(std::string &str, const std::string &from, const std::string &to, const size_t maxPosition) {
      if(from.empty())
        return;
      size_t start_pos = 0;
      size_t cur_max_pos = maxPosition;
      const size_t from_len = from.length();
      const size_t to_len = to.length();
      while((start_pos = str.find(from, start_pos)) != std::string::npos /* Loop checking, 10.07.2015, D.Mancusi */
            && (cur_max_pos==std::string::npos || start_pos<cur_max_pos)) {
        str.replace(start_pos, from_len, to);
        start_pos += to_len; // In case 'to' contains 'from', like replacing 'x' with 'yx'
        if(cur_max_pos!=std::string::npos)
          cur_max_pos += to_len - from_len;
      }
    }

    std::vector<std::string> tokenize(std::string const &str, const std::string &delimiters) {
      size_t startPos = 0, endPos;
      std::vector<std::string> tokens;
      do {
        endPos = str.find_first_of(delimiters, startPos);
        std::string token = str.substr(startPos, endPos-startPos);
        tokens.push_back(token);
        startPos = str.find_first_not_of(delimiters, endPos);
      } while(endPos!=std::string::npos); /* Loop checking, 10.07.2015, D.Mancusi */

      return tokens;
    }

    G4bool isInteger(std::string const &str) {
      const size_t pos = str.find_first_not_of("0123456789");
      return (pos==std::string::npos);
    }

    std::string expandPath(std::string const &path) {
#ifdef HAVE_WORDEXP
      wordexp_t expansionResult;
      std::string result;
      G4int err = wordexp(path.c_str(), &expansionResult, WRDE_NOCMD);
      if(err)
        result = path;
      else
        result = expansionResult.we_wordv[0];
      wordfree(&expansionResult);
      return result;
#else
      // no-op if wordexp.h is not found
      return path;
#endif
    }
  }

#endif // INCLXX_IN_GEANT4_MODE

}
