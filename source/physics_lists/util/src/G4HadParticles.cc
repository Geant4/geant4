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
// Geant4 class G4HadParticles
//
// Author V.Ivanchenko 09.05.2020
//

#include "G4HadParticles.hh"

const std::vector<G4int>& G4HadParticles::GetLightHadrons()
{
  // p, n, pi+, pi-
  static const std::vector<G4int> sLightHadrons = {
    2212, 2112, 211, -211
  };
  return sLightHadrons;
}

const std::vector<G4int>& G4HadParticles::GetHyperons()
{
  // Lambda, Sigma+, Sigma-, Xi0, Xi-, Omega-
  // (note that Sigma0 has not been included because it decays very quickly)
  static const std::vector<G4int> sHyperons = {
    3122, 3222, 3112, 3322, 3312, 3334
  };
  return sHyperons;
}

const std::vector<G4int>& G4HadParticles::GetAntiHyperons()
{
  // anti_Lambda, anti_Sigma+, anti_Sigma-, anti_Xi0, anti_Xi-, anti_Omega-
  // (note that anti_Sigma0 has not been included because it decays very quickly)
  static const std::vector<G4int> sAntiHyperons = {
    -3122, -3222, -3112, -3322, -3312, -3334
  };
  return sAntiHyperons;
}

const std::vector<G4int>& G4HadParticles::GetKaons()
{
  // K+, K-, KS, KL
  static const std::vector<G4int> sKaons = {
    321, -321, 310, 130
  };
  return sKaons;
}

const std::vector<G4int>& G4HadParticles::GetBCHadrons()
{
  // Note: etac, JPsi, SigmaC++, SigmaC+, SigmaC0, Upsilon,
  //       SigmaB+, SigmaB0, SigmaB- are not included because
  // they decay very quickly (therefore, their hadronic
  // interactions can be neglected, as for pi0 and Sigma0).
  static const std::vector<G4int> sBCHadrons = {
    // D+, D0, D-, D0bar, Ds+, Ds-
    411, 421, -411, -421, 431, -431,
    // B+, B0, B-, B0bar, Bs0, Bs0bar, Bc+, Bc-,
    521, 511, -521, -511, 531, -531, 541, -541,
    // LambdaC+, XiC+, XiC0, OmegaC0
    4122, 4232, 4132, 4332,
    // LambdaB, XiB0, XiB-, OmegaB-
    5122, 5232, 5132, 5332,
    // corresponding anti_baryons
    -4122, -4232, -4132, -4332,
    -5122, -5232, -5132, -5332
  };
  return sBCHadrons;
}

const std::vector<G4int>& G4HadParticles::GetLightIons()
{
  // d, t, He3, alpha
  static const std::vector<G4int> sLightIons = {
    1000010020, 1000010030, 1000020030, 1000020040
  };
  return sLightIons;
}

const std::vector<G4int>& G4HadParticles::GetLightAntiIons()
{
  // pbar, nbar, light anti-ions
  static const std::vector<G4int> sLightAntiIons = {
    -2212, -2112, -1000010020, -1000010030, -1000020030, -1000020040
  };
  return sLightAntiIons;
}

const std::vector<G4int>& G4HadParticles::GetHyperNuclei()
{
  // hyper_t, hyper_H4, hyper_He4, hyder_He5, 2-hyper-2n, 2-hyper_H4
  static const std::vector<G4int> sHyperNuclei = {
    1010010030, 1010010040, 1010020040, 1010020050, 1020000040, 1020010040
  };
  return sHyperNuclei;
}

const std::vector<G4int>& G4HadParticles::GetHyperAntiNuclei()
{
  // anti-hyper-nuclei
  static const std::vector<G4int> sHyperAntiNuclei = {
    -1010010030, -1010010040, -1010020040, -1010020050, -1020000040, -1020010040
  };
  return sHyperAntiNuclei;
}

const std::vector<G4int>& G4HadParticles::GetHeavyChargedParticles()
{
  //
  // charged particles for EM physics
  //
  static const std::vector<G4int> sHeavyChargedPart = {
    // Sigma+, Sigma-, Xi-, Omega-, anti_hyperons
    3222, 3112, 3312, 3334, -3222, -3112, -3312, -3334,
    // light anti_ions
    -1000010020, -1000010030, -1000020030, -1000020040,
    // tau+-
    15, -15
  };
  return sHeavyChargedPart;
}

const std::vector<G4int>& G4HadParticles::GetBCChargedHadrons()
{
  static const std::vector<G4int> sBCChargedHadrons = {
    // D+, D-, Ds+, Ds-
    411, -411, 431, -431,
    // B+, B-, Bc+, Bc-,
    521, -521, 541, -541,
    // LambdaC+, SigmaC++, SigmaC+, XiC+
    4122, 4222, 4212, 4232,
    // SigmaB+, SigmaB-, XiB-, OmegaB-
    5222, 5112, 5132, 5332,
    // anti_baryons
    -4122, -4222, -4212, -4232, -5222, -5112, -5132, -5332
  };
  return sBCChargedHadrons;
}

const std::vector<G4int>& G4HadParticles::GetChargedHyperNuclei()
{
  // hyper_t
  static const std::vector<G4int> sChargedHyperNuclei = {
    1010010030,  1010010040,  1010020040,  1010020050,  1020010040,
    -1010010030, -1010010040, -1010020040, -1010020050, -1020010040
  };
  return sChargedHyperNuclei;
}
