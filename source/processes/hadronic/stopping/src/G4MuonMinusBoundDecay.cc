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
// $Id$
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

  G4double time = -std::log(G4UniformRand()) / lambda;
  G4HadProjectile* p = const_cast<G4HadProjectile*>(&projectile);
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
  // Initialized data
  static G4double zeff[101] = { 0.,
    1.,1.98,2.95,3.89,4.8,5.72,6.61,7.49,8.32,9.12,9.95,10.69,11.48,12.22,
    12.91,13.64,14.24,14.89,15.53,16.15,16.75,17.38,18.04,18.49,
    19.06,19.59,20.1,20.66,21.12,21.61,22.02,22.43,22.84,23.24,
    23.65,24.06,24.47,24.85,25.23,25.61,25.99,26.37,26.69,27.,
    27.32,27.63,27.95,28.2,28.42,28.64,28.79,29.03,29.27,29.51,
    29.75,29.99,30.2,30.36,30.53,30.69,30.85,31.01,31.18,31.34,
    31.48,31.62,31.76,31.9,32.05,32.19,32.33,32.47,32.61,32.76,
    32.94,33.11,33.29,33.46,33.64,33.81,34.21,34.18,34.,34.1,
    34.21,34.31,34.42,34.52,34.63,34.73,34.84,34.94,35.04,35.15,
    35.25,35.36,35.46,35.57,35.67,35.78 };

  // Mu- capture data from B.B.Balashov, G.Ya.Korenman, P.A.Eramgan
  // Atomizdat, 1978. (Experimental capture velocities)
  // Data for Hydrogen from Phys. Rev. Lett. 99(2007)032002
  // Data for Helium from Phys. Rep. 354(2001)243

  const size_t ListZE = 67;
  static G4int ListZExp[ListZE] = { 1, 2,
      3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
     13, 14, 15, 16, 17, 18, 19, 20, 22, 23,
     24, 25, 26, 27, 28, 31, 32, 33, 34, 37,
     38, 39, 40, 41, 42, 45, 46, 47, 48, 49,
     50, 51, 52, 53, 55, 56, 57, 58, 59, 60,
     62, 64, 65, 67, 72, 73, 74, 80, 81, 82,
     83, 90, 92, 93};

  static G4double ListCaptureVel[ListZE] = { 0.000725, 0.000356,
     0.0057, 0.010, 0.0258, 0.0371, 0.0644,
     0.0974, 0.144, 0.250,  0.386,  0.479,
     0.700,  0.849, 1.119,  1.338,  1.40, 
     1.30,   1.98,  2.45,   2.60,   3.19,
     3.29,   3.91,  4.41,   4.96,   5.74,
     5.68,   5.53,  6.06,   5.69,   6.89,
     7.25,   7.89,  8.59,  10.40,   9.22,
    10.01,  10.00, 10.88,  10.62,  11.37,
    10.68,  10.49,  9.06,  11.20,  10.98,
    10.18,  10.71, 11.44,  13.45,  12.32,
    12.22,  12.09, 12.73,  12.95,  13.03,
    12.86,  13.13, 13.39,  12.74,  13.78,
    13.02,  13.26, 13.10,  14.00,  14.70};


  // Local variables
  G4double zeff2, xmu, a2ze, r1, r2;
  G4double lambda;

  // ==  Effective charges from Ford and Wills Nucl Phys 35(1962)295.
  // ==  Untabulated charges are interpolated.
  // ==  Mu capture lifetime (Goulard and Primakoff PRC10(1974)2034.

  G4int i = Z;
  if(i > 100) { i = 100; }

  const G4double b0a = -.03;
  const G4double b0b = -.25;
  const G4double b0c = 3.24;
  const G4double t1 = 875.e-10;
  r1 = zeff[i];
  zeff2 = r1 * r1;

  // ^-4 -> ^-5 suggested by user
  xmu = zeff2 * 2.663e-5;
  a2ze = 0.5 * A / Z;
  r2 = 1.0 - xmu;
  lambda = t1 * zeff2 * zeff2 * (r2 * r2) * (1.0 - (1.0 - xmu) * .75704) *
          (a2ze * b0a + 1.0 - (a2ze - 1.0) * b0b -
          (2 * (A - Z)  + std::fabs(a2ze - 1.) ) * b0c / G4double(A * 4) );

  // == Mu capture data are taken if exist 
  for (size_t j = 0; j < ListZE; ++j) {
    if( ListZExp[j] == i + 1) {
      lambda = ListCaptureVel[j] / microsecond;
      break;
    }
  }

  return lambda;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonMinusBoundDecay::GetMuonDecayRate(G4int Z)
{
  // Decay time on K-shell 
  // N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.

  G4double lambda = 1.0;
  if(Z > 1) { 
    G4double x = Z*fine_structure_const;
    lambda -= 2.5 * x * x; 
    if( 0.5 > lambda ) { lambda = 0.5; }
  } else {

    // Published value 0.455851 - Phys. Rev. Lett. 99(2007)032002
    lambda = 1.00151;  
  }
  return lambda * 0.445164 / microsecond; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusBoundDecay::ModelDescription(std::ostream& outFile) const
{
  outFile << "Sample probabilities of mu- nuclear capture of decay"
	  << " from K-shell orbit.\n"
	  << "  Time of projectile is changed taking into account life time"
	  << " of muonic atom.\n"
	  << "  If decay is sampled primary state become stopAndKill,"
	  << " else - isAlive.\n"
          << "Based of reviews:\n"
          << "  N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.\n"
          << "  B.B.Balashov, G.Ya.Korenman, P.A.Eramgan, Atomizdat, 1978.\n"; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

