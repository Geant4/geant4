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

#include "G4INCLGlobals.hh"
#include "G4INCLParticle.hh"

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
      while(len-startPos > lineLength) {
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
      while((start_pos = str.find(from, start_pos)) != std::string::npos
            && (cur_max_pos==std::string::npos || start_pos<cur_max_pos)) {
        str.replace(start_pos, from_len, to);
        start_pos += to_len; // In case 'to' contains 'from', like replacing 'x' with 'yx'
        if(cur_max_pos!=std::string::npos)
          cur_max_pos += to_len - from_len;
      }
    }
  }
#endif // INCLXX_IN_GEANT4_MODE

}
