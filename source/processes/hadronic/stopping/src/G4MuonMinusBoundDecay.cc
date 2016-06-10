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
// $Id: G4MuonMinusBoundDecay.cc 91836 2015-08-07 07:25:54Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class header file 
//
// File name:  G4MuonMinusBoundDecay
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 24 April 2012 on base of G4MuMinusCaptureAtRest
//
// Modified:  
// 04/23/2013  K.Genser     Fixed a constant in computation of lambda
//                          as suggested by J P Miller/Y Oksuzian;
//                          Optimized and corrected lambda calculation/lookup
// 04/30/2013  K.Genser     Improved GetMuonCaptureRate extended data and lookup 
//                          to take both Z & A into account
//                          Improved GetMuonDecayRate by using Zeff instead of Z
//                          Extracted Zeff into GetMuonZeff
//                          
//----------------------------------------------------------------------

#include "G4MuonMinusBoundDecay.hh"
#include "Randomize.hh" 
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusBoundDecay::G4MuonMinusBoundDecay()
  : G4HadronicInteraction("muMinusBoundDeacy")
{
  fMuMass = G4MuonMinus::MuonMinus()->GetPDGMass(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusBoundDecay::~G4MuonMinusBoundDecay()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadFinalState* 
G4MuonMinusBoundDecay::ApplyYourself(const G4HadProjectile& projectile, 
                                     G4Nucleus& targetNucleus)
{
  result.Clear();
  G4int Z = targetNucleus.GetZ_asInt(); 
  G4int A = targetNucleus.GetA_asInt(); 

  // Decide on Decay or Capture, and doit.
  G4double lambdac  = GetMuonCaptureRate(Z, A);
  G4double lambdad  = GetMuonDecayRate(Z);
  G4double lambda   = lambdac + lambdad;

  // ===  sample capture  time and change time of projectile
  // ===  this is needed for the case when bound decay is not happen
  // ===  but muon is capruted by the nucleus with some delay

  G4HadProjectile* p = const_cast<G4HadProjectile*>(&projectile);
  G4double time = p->GetGlobalTime() - G4Log(G4UniformRand())/lambda;
  p->SetGlobalTime(time);
    
  //G4cout << "lambda= " << lambda << " lambdac= " << lambdac 
  //<< " t= " << time << G4endl;
 
  // cascade
  if( G4UniformRand()*lambda < lambdac) {
    result.SetStatusChange(isAlive);

  } else {

    // Simulation on Decay of mu- on a K-shell of the muonic atom
    result.SetStatusChange(stopAndKill);
    G4double xmax = 1 + electron_mass_c2*electron_mass_c2/(fMuMass*fMuMass);
    G4double xmin = 2.0*electron_mass_c2/fMuMass;
    G4double KEnergy = projectile.GetBoundEnergy();

    /*
      G4cout << "G4MuonMinusBoundDecay::ApplyYourself" 
      << " XMAX= " << xmax << " Ebound= " << KEnergy<< G4endl;
    */
    G4double pmu = std::sqrt(KEnergy*(KEnergy + 2.0*fMuMass));
    G4double emu = KEnergy + fMuMass;
    G4ThreeVector dir = G4RandomDirection();
    G4LorentzVector MU(pmu*dir, emu);
    G4ThreeVector bst = MU.boostVector();

    G4double Eelect, Pelect, x, ecm;
    G4LorentzVector EL, NN;
    // Calculate electron energy
    // these do/while loops are safe
    do {
      do {
        x = xmin + (xmax-xmin)*G4UniformRand();
      } while (G4UniformRand() > (3.0 - 2.0*x)*x*x );
      Eelect = x*fMuMass*0.5;
      Pelect = 0.0;
      if(Eelect > electron_mass_c2) { 
        Pelect = std::sqrt(Eelect*Eelect - electron_mass_c2*electron_mass_c2);
      } else {
        Pelect = 0.0;
        Eelect = electron_mass_c2;
      }
      dir = G4RandomDirection();
      EL = G4LorentzVector(Pelect*dir,Eelect);
      EL.boost(bst);
      Eelect = EL.e() - electron_mass_c2 - 2.0*KEnergy;
      //
      // Calculate rest frame parameters of 2 neutrinos
      //
      NN = MU - EL;
      ecm = NN.mag2();
      // Loop checking, 06-Aug-2015, Vladimir Ivanchenko
    } while (Eelect < 0.0 || ecm < 0.0);

    //
    // Create electron
    //
    G4DynamicParticle* dp = new G4DynamicParticle(G4Electron::Electron(),
                                                  EL.vect().unit(),
                                                  Eelect);

    AddNewParticle(dp, time);
    //
    // Create Neutrinos
    //
    ecm = 0.5*std::sqrt(ecm);
    bst = NN.boostVector();
    G4ThreeVector p1 = ecm * G4RandomDirection();
    G4LorentzVector N1 = G4LorentzVector(p1,ecm);
    N1.boost(bst);
    dp = new G4DynamicParticle(G4AntiNeutrinoE::AntiNeutrinoE(), N1);
    AddNewParticle(dp, time);
    NN -= N1;
    dp = new G4DynamicParticle(G4NeutrinoMu::NeutrinoMu(), NN);
    AddNewParticle(dp, time);
  }
  return &result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonMinusBoundDecay::GetMuonCaptureRate(G4int Z, G4int A)
{

  // Initialize data

  // Mu- capture data from 
  // T. Suzuki, D. F. Measday, J.P. Roalsvig Phys.Rev. C35 (1987) 2212
  // weighted average of the two most precise measurements

  // Data for Hydrogen from Phys. Rev. Lett. 99(2007)032002
  // Data for Helium from D.F. Measday Phys. Rep. 354(2001)243

  struct capRate {
    G4int        Z;
    G4int        A;
    G4double cRate;
    G4double cRErr;
  };

  // this struct has to be sorted by Z when initialized as we exit the
  // loop once Z is above the stored value; cRErr are not used now but
  // are included for completeness and future use

  const capRate capRates [] = {
    {  1,   1,  0.000725, 0.000017 },
    {  2,   3,  0.002149, 0.00017 }, 
    {  2,   4,  0.000356, 0.000026 },
    {  3,   6,  0.004647, 0.00012 },
    {  3,   7,  0.002229, 0.00012 },
    {  4,   9,  0.006107, 0.00019 },
    {  5,  10,  0.02757 , 0.00063 },
    {  5,  11,  0.02188 , 0.00064 },
    {  6,  12,  0.03807 , 0.00031 },
    {  6,  13,  0.03474 , 0.00034 },
    {  7,  14,  0.06885 , 0.00057 },
    {  8,  16,  0.10242 , 0.00059 },
    {  8,  18,  0.0880  , 0.0015  },
    {  9,  19,  0.22905 , 0.00099 },
    { 10,  20,  0.2288  , 0.0045 },
    { 11,  23,  0.3773  , 0.0014 },
    { 12,  24,  0.4823  , 0.0013 },
    { 13,  27,  0.6985  , 0.0012 },
    { 14,  28,  0.8656  , 0.0015 },
    { 15,  31,  1.1681  , 0.0026 },
    { 16,  32,  1.3510  , 0.0029 },
    { 17,  35,  1.800   , 0.050 },
    { 17,  37,  1.250   , 0.050 },
    { 18,  40,  1.2727  , 0.0650 },
    { 19,  39,  1.8492  , 0.0050 },
    { 20,  40,  2.5359  , 0.0070 },
    { 21,  45,  2.711   , 0.025 },
    { 22,  48,  2.5908  , 0.0115 },
    { 23,  51,  3.073   , 0.022 },
    { 24,  50,  3.825   , 0.050 },
    { 24,  52,  3.465   , 0.026 },
    { 24,  53,  3.297   , 0.045 },
    { 24,  54,  3.057   , 0.042 },
    { 25,  55,  3.900   , 0.030 },
    { 26,  56,  4.408   , 0.022 },
    { 27,  59,  4.945   , 0.025 },
    { 28,  58,  6.11    , 0.10 },
    { 28,  60,  5.56    , 0.10 },
    { 28,  62,  4.72    , 0.10 },
    { 29,  63,  5.691   , 0.030 },
    { 30,  66,  5.806   , 0.031 },
    { 31,  69,  5.700   , 0.060 },
    { 32,  72,  5.561   , 0.031 },
    { 33,  75,  6.094   , 0.037 },
    { 34,  80,  5.687   , 0.030 },
    { 35,  79,  7.223   , 0.28 },
    { 35,  81,  7.547   , 0.48 },
    { 37,  85,  6.89    , 0.14 },
    { 38,  88,  6.93    , 0.12 },
    { 39,  89,  7.89    , 0.11 },
    { 40,  91,  8.620   , 0.053 },
    { 41,  93, 10.38    , 0.11 },
    { 42,  96,  9.298   , 0.063 },
    { 45, 103, 10.010   , 0.045 },
    { 46, 106, 10.000   , 0.070 },
    { 47, 107, 10.869   , 0.095 },
    { 48, 112, 10.624   , 0.094 },
    { 49, 115, 11.38    , 0.11 },
    { 50, 119, 10.60    , 0.11 },
    { 51, 121, 10.40    , 0.12 },
    { 52, 128,  9.174   , 0.074 },
    { 53, 127, 11.276   , 0.098 },
    { 55, 133, 10.98    , 0.25 },
    { 56, 138, 10.112   , 0.085 },
    { 57, 139, 10.71    , 0.10 },
    { 58, 140, 11.501   , 0.087 },
    { 59, 141, 13.45    , 0.13 },
    { 60, 144, 12.35    , 0.13 },
    { 62, 150, 12.22    , 0.17 },
    { 64, 157, 12.00    , 0.13 },
    { 65, 159, 12.73    , 0.13 },
    { 66, 163, 12.29    , 0.18 },
    { 67, 165, 12.95    , 0.13 },
    { 68, 167, 13.04    , 0.27 },
    { 72, 178, 13.03    , 0.21 },
    { 73, 181, 12.86    , 0.13 },
    { 74, 184, 12.76    , 0.16 },
    { 79, 197, 13.35    , 0.10 },
    { 80, 201, 12.74    , 0.18 },
    { 81, 205, 13.85    , 0.17 },
    { 82, 207, 13.295   , 0.071 },
    { 83, 209, 13.238   , 0.065 },
    { 90, 232, 12.555   , 0.049 },
    { 92, 238, 12.592   , 0.035 },
    { 92, 233, 14.27    , 0.15 },
    { 92, 235, 13.470   , 0.085 },
    { 92, 236, 13.90    , 0.40 },
    { 93, 237, 13.58    , 0.18 },
    { 94, 239, 13.90    , 0.20 },
    { 94, 242, 12.86    , 0.19 }
  };

  G4double lambda = -1.;

  size_t nCapRates = sizeof(capRates)/sizeof(capRates[0]);
  for (size_t j = 0; j < nCapRates; ++j) {
    if( capRates[j].Z == Z && capRates[j].A == A ) {
      lambda = capRates[j].cRate / microsecond;
      break;
    }
    // make sure the data is sorted for the next statement to work correctly
    if (capRates[j].Z > Z) {break;}
  }

  if (lambda < 0.) {

    // ==  Mu capture lifetime (Goulard and Primakoff PRC10(1974)2034.

    const G4double b0a = -0.03;
    const G4double b0b = -0.25;
    const G4double b0c =  3.24;
    const G4double t1 = 875.e-9; // -10-> -9  suggested by user
    G4double r1 = GetMuonZeff(Z);
    G4double zeff2 = r1 * r1;

    // ^-4 -> ^-5 suggested by user
    G4double xmu = zeff2 * 2.663e-5;
    G4double a2ze = 0.5 *G4double(A) / G4double(Z);
    G4double r2 = 1.0 - xmu;
    lambda = t1 * zeff2 * zeff2 * (r2 * r2) * (1.0 - (1.0 - xmu) * .75704) *
      (a2ze * b0a + 1.0 - (a2ze - 1.0) * b0b -
       G4double(2 * (A - Z)  + std::abs(a2ze - 1.) ) * b0c / G4double(A * 4) );

  }

  return lambda;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4MuonMinusBoundDecay::GetMuonZeff(G4int Z)
{

  // ==  Effective charges from 
  // "Total Nuclear Capture Rates for Negative Muons"
  // T. Suzuki, D. F. Measday, J.P. Roalsvig Phys.Rev. C35 (1987) 2212
  // and if not present from
  // Ford and Wills Nucl Phys 35(1962)295 or interpolated


  const size_t maxZ = 100;
  const G4double zeff[maxZ+1] = 
    {  0.,
       1.00, 1.98, 2.94, 3.89, 4.81, 5.72, 6.61, 7.49, 8.32, 9.14,
       9.95,10.69,11.48,12.22,12.90,13.64,14.24,14.89,15.53,16.15,
       16.77,17.38,18.04,18.49,19.06,19.59,20.13,20.66,21.12,21.61,
       22.02,22.43,22.84,23.24,23.65,24.06,24.47,24.85,25.23,25.61,
       25.99,26.37,26.69,27.00,27.32,27.63,27.95,28.20,28.42,28.64,
       28.79,29.03,29.27,29.51,29.75,29.99,30.22,30.36,30.53,30.69,
       30.85,31.01,31.18,31.34,31.48,31.62,31.76,31.90,32.05,32.19,
       32.33,32.47,32.61,32.76,32.94,33.11,33.29,33.46,33.64,33.81,
       34.21,34.18,34.00,34.10,34.21,34.31,34.42,34.52,34.63,34.73,
       34.84,34.94,35.05,35.16,35.25,35.36,35.46,35.57,35.67,35.78 };

  if (Z<0) {Z=0;}
  if (Z>G4int(maxZ)) {Z=maxZ;}

  return zeff[Z];

}


G4double G4MuonMinusBoundDecay::GetMuonDecayRate(G4int Z)
{
  // Decay time on K-shell 
  // N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.

  // this is the "small Z" approximation formula (2.9)
  // Lambda(bound)/Lambda(free) = 1-beta(Z*alpha)**2 with beta~=2.5
  // we assume that Z is Zeff

  // PDG 2012 muon lifetime value is 2.1969811(22) 10e-6s
  // which when inverted gives       0.45517005    10e+6/s

  struct decRate {
    G4int         Z;
    G4double  dRate;
    G4double  dRErr;
  };

  // this struct has to be sorted by Z when initialized as we exit the
  // loop once Z is above the stored value

  const decRate decRates [] = {
    {  1,  0.4558514, 0.0000151 }
  };

  G4double lambda = -1.;

  // size_t nDecRates = sizeof(decRates)/sizeof(decRates[0]);
  // for (size_t j = 0; j < nDecRates; ++j) {
  //   if( decRates[j].Z == Z ) {
  //     lambda = decRates[j].dRate / microsecond;
  //     break;
  //   }
  //   // make sure the data is sorted for the next statement to work
  //   if (decRates[j].Z > Z) {break;}
  // }

  // we'll use the above code once we have more data
  // since we only have one value we just assign it
  if (Z == 1) {lambda =  decRates[0].dRate/microsecond;}

  if (lambda < 0.) {
    const G4double freeMuonDecayRate =  0.45517005 / microsecond;
    lambda = 1.0;
    G4double x = GetMuonZeff(Z)*fine_structure_const;
    lambda -= 2.5 * x * x;
    lambda *= freeMuonDecayRate;
  }

  return lambda;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusBoundDecay::ModelDescription(std::ostream& outFile) const
{
  outFile << "Sample probabilities of mu- nuclear capture of decay"
          << " from K-shell orbit.\n"
          << "   Time of projectile is changed taking into account life time"
          << " of muonic atom.\n"
          << "   If decay is sampled primary state become stopAndKill,"
          << " else - isAlive.\n"
          << "  Based of reviews:\n"
          << "  N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.\n"
          << "  T. Suzuki, D. F. Measday, J.P. Roalsvig Phys.Rev. C35 (1987) 2212\n"; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

