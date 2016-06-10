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

/** \file G4INCLHornerFormEvaluator.hh
 * \brief Template-metaprogramming-based evaluator for polynomials in Horner form
 *
 * \date 3rd March 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLHORNERFORMEVALUATOR_HH
#define G4INCLHORNERFORMEVALUATOR_HH

namespace G4INCL {

  template<G4int N>
    class HornerCoefficients {
      protected:
        G4double a[N];
      public:
        G4double &operator[](G4int i) { return a[i]; }
        const G4double &operator[](G4int i) const { return a[i]; }
    };

  struct HornerC1 : public HornerCoefficients<1> {
    HornerC1(
             const G4double a0
            ) {
      a[0] = a0;
    }
  };

  struct HornerC2 : public HornerCoefficients<2> {
    HornerC2(
             const G4double a0,
             const G4double a1
            ) {
      a[0] = a0;
      a[1] = a1;
    }
  };

  struct HornerC3 : public HornerCoefficients<3> {
    HornerC3(
             const G4double a0,
             const G4double a1,
             const G4double a2
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
    }
  };

  struct HornerC4 : public HornerCoefficients<4> {
    HornerC4(
             const G4double a0,
             const G4double a1,
             const G4double a2,
             const G4double a3
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
    }
  };

  struct HornerC5 : public HornerCoefficients<5> {
    HornerC5(
             const G4double a0,
             const G4double a1,
             const G4double a2,
             const G4double a3,
             const G4double a4
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
    }
  };

  struct HornerC6 : public HornerCoefficients<6> {
    HornerC6(
             const G4double a0,
             const G4double a1,
             const G4double a2,
             const G4double a3,
             const G4double a4,
             const G4double a5
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
      a[5] = a5;
    }
  };

  struct HornerC7 : public HornerCoefficients<7> {
    HornerC7(
             const G4double a0,
             const G4double a1,
             const G4double a2,
             const G4double a3,
             const G4double a4,
             const G4double a5,
             const G4double a6
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
      a[5] = a5;
      a[6] = a6;
    }
  };

  struct HornerC8 : public HornerCoefficients<8> {
    HornerC8(
             const G4double a0,
             const G4double a1,
             const G4double a2,
             const G4double a3,
             const G4double a4,
             const G4double a5,
             const G4double a6,
             const G4double a7
            ) {
      a[0] = a0;
      a[1] = a1;
      a[2] = a2;
      a[3] = a3;
      a[4] = a4;
      a[5] = a5;
      a[6] = a6;
      a[7] = a7;
    }
  };

  template<G4int M>
    struct HornerEvaluator {
      template<G4int N>
        static G4double eval(const G4double x, HornerCoefficients<N> const &coeffs) {
          return coeffs[N-M] + x * HornerEvaluator<M-1>::eval(x, coeffs);
        }
    };

  template<>
    struct HornerEvaluator<1> {
      template<G4int N>
        static G4double eval(const G4double, HornerCoefficients<N> const &coeffs) {
          return coeffs[N-1];
        }
    };

}

#endif
