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

/*
 * G4INCLRandom.cc
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
// #include <cassert>

namespace G4INCL {

  namespace Random {

    namespace {

      G4ThreadLocal IRandomGenerator* theGenerator;

    }

    void setGenerator(G4INCL::IRandomGenerator *aGenerator) {
      if(isInitialized()) {
        INCL_ERROR("INCL random number generator already initialized." << std::endl);
      } else {
        theGenerator = aGenerator;
      }
    }

    void setSeeds(const SeedVector &sv) {
      theGenerator->setSeeds(sv);
    }

    SeedVector getSeeds() {
      return theGenerator->getSeeds();
    }

    G4double shoot() {
      return theGenerator->flat();
    }

    G4double shoot0() {
      G4double r;
      while( (r=shoot()) <= 0. )
        ;
      return r;
    }

    G4double shoot1() {
      G4double r;
      while( (r=shoot()) >= 1. )
        ;
      return r;
    }

    template<typename T>
      T shootInteger(T n) {
        return static_cast<T>(shoot1() * n);
      }

    G4double gauss(G4double sigma) {
      // generate a Gaussian random number with standard deviation sigma
      // uses the flat() and flat0() methods
      static G4ThreadLocal G4bool generated = false;
      static G4ThreadLocal G4double u, v;

      if( !generated )
      {
        u = shoot0();
        v = Math::twoPi*shoot();
        generated = true;
        return sigma*std::sqrt(-2*std::log(u))*std::cos(v);
      }
      else
      {
        generated = false;
        return sigma*std::sqrt(-2*std::log(u))*std::sin(v);
      }
    }

    ThreeVector normVector(G4double norm) {

      const G4double ctheta = (1.-2.*shoot());
      const G4double stheta = std::sqrt(1.-ctheta*ctheta);
      const G4double phi = Math::twoPi*shoot();
      return ThreeVector(
                         norm * stheta * std::cos(phi),
                         norm * stheta * std::sin(phi),
                         norm * ctheta);

    }

    ThreeVector sphereVector(G4double rmax) {
      return normVector( rmax*Math::pow13(shoot0()) );
    }

    ThreeVector gaussVector(G4double sigma) {
      const G4double sigmax = sigma * Math::oneOverSqrtThree;
      return ThreeVector(gauss(sigmax), gauss(sigmax), gauss(sigmax));
    }

    std::pair<G4double,G4double> correlatedGaussian(const G4double corrCoeff, const G4double x0, const G4double sigma) {
// assert(corrCoeff<=1. && corrCoeff>=-1.);
      G4double factor = 1.-corrCoeff*corrCoeff;
      if(factor<=0.)
        factor=0.;
      const G4double x = gauss(sigma) + x0;
      const G4double y = corrCoeff * x + gauss(sigma*std::sqrt(factor)) + x0;
      return std::make_pair(x, y);
    }

    std::pair<G4double,G4double> correlatedUniform(const G4double corrCoeff) {
      std::pair<G4double,G4double> gaussians = correlatedGaussian(corrCoeff);
      return std::make_pair(Math::gaussianCDF(gaussians.first), Math::gaussianCDF(gaussians.second));
    }

    void deleteGenerator() {
      delete theGenerator;
      theGenerator = 0;
    }

    G4bool isInitialized() {
      if(theGenerator == 0) return false;
      return true;
    }

  }

}
