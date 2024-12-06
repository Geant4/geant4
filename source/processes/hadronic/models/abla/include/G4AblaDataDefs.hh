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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, UDC (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (initial translation of ablav3p)
// Aleksandra Kelic, GSI (ABLA07 code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//

#pragma once

#include "globals.hh"
#include <cmath>
#include <vector>


constexpr const G4int nrows = 180;
constexpr const G4int zcols = 122;

constexpr const G4int lpcols = 13;
constexpr const G4int lprows = 154;

constexpr const G4int nrowsbeta = 251;
constexpr const G4int zcolsbeta = 137;

constexpr const G4int indexpart = 300;

// Data structures needed by ABLA evaporation code

class G4Mexp {

public:
  G4Mexp(){};

  virtual ~G4Mexp() = default;

  G4double massexp[lprows][lpcols] = {{0.}};
  G4double bind[lprows][lpcols] = {{0.}};
  G4int mexpiop[lprows][lpcols] = {{0}};
};

class G4Ec2sub {
public:
  G4Ec2sub(){};

  virtual ~G4Ec2sub() = default;

  G4double ecnz[nrows][zcols] = {{0.}};
};

class G4Ald {
public:
  G4Ald() : av(0.0), as(0.0), ak(0.0), optafan(0.0){};

  virtual ~G4Ald() = default;

  G4double av, as, ak, optafan = 0.;
};

/**
 * Shell corrections and deformations.
 **/

class G4Ecld {

public:
  G4Ecld(){};
  virtual ~G4Ecld() = default;

  /**
   * Ground state shell correction frldm for a spherical ground state.
   */
  G4double ecgnz[nrows][zcols] = {{0.}};

  /**
   * Shell correction for the saddle point (now: == 0).
   */
  G4double ecfnz[nrows][zcols] = {{0.}};

  /**
   * Difference between deformed ground state and ldm value.
   */
  G4double vgsld[nrows][zcols] = {{0.}};

  /**
   * Alpha ground state deformation (this is not beta2!)
   * beta2 = std::sqrt(5/(4pi)) * alpha
   */
  G4double alpha[nrows][zcols] = {{0.}};

  /**
   * RMS function for lcp emission barriers
   */
  G4double rms[nrows][zcols] = {{0.}};

  /**
   * Beta2 deformations
   */
  G4double beta2[nrowsbeta][zcolsbeta] = {{0.}};

  /**
   * Beta4 deformations
   */
  G4double beta4[nrowsbeta][zcolsbeta] = {{0.}};
};

class G4Fiss {
  /**
   * Options and parameters for fission channel.
   */

public:
  G4Fiss()
      : bet(0.0), bethyp(0.0), ifis(0.0), ucr(0.0), dcr(0.0), optshp(0), optxfis(0),
        optct(0), optcol(0), at(0), zt(0){};

  virtual ~G4Fiss() = default;

  G4double bet, bethyp, ifis, ucr, dcr;
  G4int optshp, optxfis, optct, optcol, at, zt;
};

/**
 * Fission barriers.
 */

class G4Fb {

public:
  G4Fb(){};
  virtual ~G4Fb() = default;

  G4double efa[nrows][zcols] = {{0.}};
};

/**
 * Options
 */

class G4Opt {

public:
  G4Opt() : optemd(0), optcha(0), optshpimf(0), optimfallowed(0), nblan0(0){};
  
  virtual ~G4Opt() = default;

  G4int optemd, optcha, optshpimf, optimfallowed, nblan0;
};

class G4VarNtp {
public:
  G4VarNtp() { clear(); };

  virtual ~G4VarNtp() = default;

  void clear() {
    ntrack = 0;
    kfis = 0;
    itypcasc.clear();
    avv.clear();
    zvv.clear();
    svv.clear();
    enerj.clear();
    pxlab.clear();
    pylab.clear();
    pzlab.clear();
  }

  /**
   * Fission 1/0=Y/N.
   */
  G4int kfis;

  /**
   * Excit energy at fis.
   */
  G4double estfis = 0.;

  /**
   * Z of fiss nucleus.
   */
  G4int izfis = 0;

  /**
   * A of fiss nucleus.
   */
  G4int iafis = 0;

  /**
   * Number of particles.
   */
  G4int ntrack;

  /**
   * Does this nucleus require Fermi break-up treatment? Only
   * applicable when used together with Geant4.
   * true = do fermi break-up (and skip ABLA part)
   * false = use ABLA
   */
  G4bool needsFermiBreakup = false;

  /**
   * emitted in cascade (0) or evaporation (1).
   */
  std::vector<G4int> itypcasc;

  /**
   * A (-1 for pions).
   */
  std::vector<G4int> avv;

  /**
   * Z
   */
  std::vector<G4int> zvv;

  /**
   * S (-1 for lambda_0).
   */
  std::vector<G4int> svv;

  /**
   * Kinetic energy.
   */
  std::vector<G4double> enerj;

  /**
   * Momentum.
   */
  std::vector<G4double> plab;
  std::vector<G4double> pxlab;
  std::vector<G4double> pylab;
  std::vector<G4double> pzlab;

  /**
   * Theta angle.
   */
  std::vector<G4double> tetlab;

  /**
   * Phi angle.
   */
  std::vector<G4double> philab;

};
