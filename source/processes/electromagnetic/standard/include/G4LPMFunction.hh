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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name: G4LPMFunction
//
// Author: Mihaly Novak
//
// Creation date: 08 October 2025
//
// Modifications:
//
// Class Description:
//
// The `G(s)` and `\phi(s)` Landau Pomeranchuk Migdal (LPM) suppression
// functions, utilised both in some bremsstrahlung and pair-production models
// when computing the actual LPM suppression effect, are factored in this util.
// The functions are pre-computed and interpolated over the `s \in [0,2)`
// intervall while approximated at `s > 2` values when they converge to unity.
//
// `G(s)` and `\phi(s)` are the spin flip and no flip suppression functions of
// Migdal (Migdal PR,1956) including some slowly converging series to which
// approximate expressions were derived by Stanev (Stanev at al. PRD,1982).
// (See e.g. `G4eBremsstrahlungRelModel::ComputeLPMGsPhis` at Geant4-v.11.3.0
//  on the details of the computation.)
//

// -------------------------------------------------------------------
//

#ifndef G4LPMFunction_h
#define G4LPMFunction_h 1

#include "G4Types.hh"

namespace G4LPMFunction {
  // Precomputed G(s) and Phi(s) LPM functions.
  // Grid: s \in [0, 2.0], ds = 0.05 -> N = 41 points, interleaved (G(s),Phi(s))
  inline constexpr G4double kFuncLPM[] = {
    0.0000E+00, 0.0000E+00,  6.9163E-02, 2.5747E-01,  2.0597E-01, 4.4573E-01,
    3.5098E-01, 5.8373E-01,  4.8095E-01, 6.8530E-01,  5.8926E-01, 7.6040E-01,
    6.7626E-01, 8.1626E-01,  7.4479E-01, 8.5805E-01,  7.9826E-01, 8.8952E-01,
    8.4003E-01, 9.1338E-01,  8.7258E-01, 9.3159E-01,  8.9794E-01, 9.4558E-01,
    9.1776E-01, 9.5640E-01,  9.3332E-01, 9.6483E-01,  9.4560E-01, 9.7143E-01,
    9.5535E-01, 9.7664E-01,  9.6313E-01, 9.8078E-01,  9.6939E-01, 9.8408E-01,
    9.7444E-01, 9.8673E-01,  9.7855E-01, 9.8888E-01,  9.8191E-01, 9.9062E-01,
    9.8467E-01, 9.9204E-01,  9.8695E-01, 9.9321E-01,  9.8884E-01, 9.9417E-01,
    9.9042E-01, 9.9497E-01,  9.9174E-01, 9.9564E-01,  9.9285E-01, 9.9619E-01,
    9.9379E-01, 9.9666E-01,  9.9458E-01, 9.9706E-01,  9.9526E-01, 9.9739E-01,
    9.9583E-01, 9.9768E-01,  9.9632E-01, 9.9794E-01,  9.9674E-01, 9.9818E-01,
    9.9710E-01, 9.9839E-01,  9.9741E-01, 9.9857E-01,  9.9767E-01, 9.9873E-01,
    9.9790E-01, 9.9887E-01,  9.9809E-01, 9.9898E-01,  9.9826E-01, 9.9909E-01,
    9.9840E-01, 9.9918E-01,  9.9856E-01, 9.9926E-01
  };

  // Obtain the `G(s)` and `\phi(s)` LPM suppression functions at any `s >= 0`.
  inline void GetLPMFunctions(G4double& lpmFuncG, G4double& lpmFuncPhi, G4double sVar) {
    // sanity check (s should be >= 0)
    if (sVar < 0.0) {
      lpmFuncG = 0.0;
      lpmFuncPhi = 0.0;
      return;
    }
    // case of `s in [0, 2)` use the precomputed functions and interpolate
    const G4double lpmSLimit =  2.0; // max_s:=2
    const G4double lpmISDelt = 20.0; // deta_s:=0.05, 1/delta_s=20
    if (sVar < lpmSLimit) {
      G4double  val = sVar*lpmISDelt;
      G4int    ilow = static_cast<G4int>(val);
      val       -= ilow;
      ilow      *= 2;
      lpmFuncG   = (kFuncLPM[ilow+2] - kFuncLPM[ilow]  )*val + kFuncLPM[ilow];
      lpmFuncPhi = (kFuncLPM[ilow+3] - kFuncLPM[ilow+1])*val + kFuncLPM[ilow+1];
      return;
    }
    // asymptotic case: G(s), Phi(s) goes to 1.0
    G4double ss = 1.0/(sVar*sVar);
    ss *= ss;
    lpmFuncG   = 1.0 - 0.0230655*ss;
    lpmFuncPhi = 1.0 - 0.01190476*ss;
  }
};

#endif //  G4LPMFunction_h
