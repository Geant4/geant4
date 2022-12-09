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
// GEANT4 Class file
//
//
// File name:     G4MuonToMuonPairProductionModel
//
// Author:        Siddharth Yajaman on the base of Vladimir Ivantchenko code
//
// Creation date: 12.07.2022
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuonToMuonPairProductionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// static members
//
static const G4int nzdat = 5;
static const G4int zdat[5] = {1, 4, 13, 29, 92};

static const G4double xgi[] =
{ 0.0198550717512320, 0.1016667612931865, 0.2372337950418355, 0.4082826787521750,
  0.5917173212478250, 0.7627662049581645, 0.8983332387068135, 0.9801449282487680 };

static const G4double wgi[] =
{ 0.0506142681451880, 0.1111905172266870, 0.1568533229389435, 0.1813418916891810,
  0.1813418916891810, 0.1568533229389435, 0.1111905172266870, 0.0506142681451880 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuonToMuonPairProductionModel::G4MuonToMuonPairProductionModel(
                                 const G4ParticleDefinition* p,
                                 const G4String& nam)
  : G4MuPairProductionModel(p, nam)
{
  theMuonMinus = G4MuonMinus::MuonMinus();
  theMuonPlus = G4MuonPlus::MuonPlus();
  muonMass = theMuonPlus->GetPDGMass();
  minPairEnergy = 2.*muonMass;
  mueRatio = muonMass/CLHEP::electron_mass_c2;
  factorForCross = 2./(3*CLHEP::pi)*
    pow(CLHEP::fine_structure_const*CLHEP::classic_electr_radius/mueRatio, 2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuonToMuonPairProductionModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy)
// Calculates the  differential (D) microscopic cross section
// using the cross section formula of Kelner, Kokoulin and Petrukhin (1999)
// Code written by Siddharth Yajaman (12/07/2022)
{
  if (pairEnergy <= minPairEnergy)
    return 0.0;

  G4double totalEnergy = tkin + particleMass;
  G4double residEnergy = totalEnergy - pairEnergy;

  if (residEnergy <= muonMass) { return 0.0; }

  G4double a0 = 1.0 / (totalEnergy * residEnergy);
  G4double rhomax = 1.0 - 2*muonMass/pairEnergy;
  G4double tmnexp = 1. - rhomax;

  if(tmnexp >= 1.0) { return 0.0; }

  G4double tmn = G4Log(tmnexp);

  G4double z2 = Z*Z;
  G4double beta = 0.5*pairEnergy*pairEnergy*a0;
  G4double xi0 = 0.5*beta;

  // Gaussian integration in ln(1-ro) ( with 8 points)
  G4double rho[8];
  G4double rho2[8];
  G4double xi[8];
  G4double xi1[8];
  G4double xii[8];

  for (G4int i = 0; i < 8; ++i)
  {
    rho[i] = G4Exp(tmn*xgi[i]) - 1.0; // rho = -asymmetry
    rho2[i] = rho[i] * rho[i];
    xi[i] = xi0*(1.0-rho2[i]);
    xi1[i] = 1.0 + xi[i];
    xii[i] = 1.0 / xi[i];
  }

  G4double ximax = xi0*(1. - rhomax*rhomax);

  G4double Y = 10 * sqrt(particleMass/totalEnergy);
  G4double U[8];

  for (G4int i = 0; i < 8; ++i)
  {
    U[i] = U_func(Z, rho2[i], xi[i], Y, pairEnergy);
  }

  G4double UMax = U_func(Z, rhomax*rhomax, ximax, Y, pairEnergy);

  G4double sum = 0.0;
  for (G4int i = 0; i < 8; ++i)
  {
    G4double X = 1 + U[i] - UMax;
    G4double lnX = G4Log(X);
    G4double phi = ((2 + rho2[i])*(1 + beta) + xi[i]*(3 + rho2[i]))*
                    G4Log(1 + xii[i]) - 1 - 3*rho2[i] + beta*(1 - 2*rho2[i])
                    + ((1 + rho2[i])*(1 + 1.5*beta) - xii[i]*(1 + 2*beta)
                    *(1 - rho2[i]))*G4Log(xi1[i]);
    sum += wgi[i]*(1.0 + rho[i])*phi*lnX;
  }

  return -tmn*sum*factorForCross*z2*residEnergy/(totalEnergy*pairEnergy);

}

G4double G4MuonToMuonPairProductionModel::U_func(G4double ZZ, G4double rho2, 
                                                 G4double xi, G4double Y,
                                                 G4double pairEnergy,
                                                 const G4double B)
{
  G4int Z = G4lrint(ZZ);
  G4double A27 = nist->GetA27(Z);
  G4double Z13 = nist->GetZ13(Z);
  static const G4double sqe = std::sqrt(G4Exp(1.0));
  G4double res = (0.65 * B / (A27*Z13) * mueRatio)/
    (1 + (2*sqe*muonMass*muonMass*(B/Z13)*(1 + xi)*(1 + Y))
     /(CLHEP::electron_mass_c2*pairEnergy*(1 - rho2)));
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuonToMuonPairProductionModel::SampleSecondaries(
                              std::vector<G4DynamicParticle*>* vdp, 
                              const G4MaterialCutsCouple* couple,
                              const G4DynamicParticle* aDynamicParticle,
                              G4double tmin,
                              G4double tmax)
{
  G4double kinEnergy = aDynamicParticle->GetKineticEnergy();
  //G4cout << "--- G4MuonToMuonPairProductionModel::SampleSecondaries E(MeV)= " 
  //         << kinEnergy << "  " 
  //         << aDynamicParticle->GetDefinition()->GetParticleName() << G4endl;
  G4double totalMomentum = std::sqrt(kinEnergy*(kinEnergy + 2.0*muonMass));

  G4ThreeVector partDirection = aDynamicParticle->GetMomentumDirection();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple,particle,kinEnergy);

  // define interval of energy transfer
  G4double maxPairEnergy = MaxSecondaryEnergyForElement(kinEnergy, 
                                                        anElement->GetZ());
  G4double maxEnergy = std::min(tmax, maxPairEnergy);
  G4double minEnergy = std::max(tmin, minPairEnergy);

  if(minEnergy >= maxEnergy) { return; }
  //G4cout << "emin= " << minEnergy << " emax= " << maxEnergy 
  // << " minPair= " << minPairEnergy << " maxpair= " << maxPairEnergy 
  //    << " ymin= " << ymin << " dy= " << dy << G4endl;

  G4double coeff = G4Log(minPairEnergy/kinEnergy)/ymin;

  // compute limits 
  G4double yymin = G4Log(minEnergy/kinEnergy)/coeff;
  G4double yymax = G4Log(maxEnergy/kinEnergy)/coeff;
 
  //G4cout << "yymin= " << yymin << "  yymax= " << yymax << G4endl;

  // units should not be used, bacause table was built without
  G4double logTkin = G4Log(kinEnergy/CLHEP::MeV);

  // sample mu-mu+ pair energy first

  // select sample table via Z
  G4int iz1(0), iz2(0);
  for(G4int iz=0; iz<nzdat; ++iz) { 
    if(currentZ == zdat[iz]) {
      iz1 = iz2 = currentZ; 
      break;
    } else if(currentZ < zdat[iz]) {
      iz2 = zdat[iz];
      if(iz > 0) { iz1 = zdat[iz-1]; }
      else { iz1 = iz2; }
      break;
    } 
  }
  if(0 == iz1) { iz1 = iz2 = zdat[nzdat-1]; }

  G4double pairEnergy = 0.0;
  G4int count = 0;
  //G4cout << "start loop Z1= " << iz1 << " Z2= " << iz2 << G4endl;
  do {
    ++count;
    // sampling using only one random number
    G4double rand = G4UniformRand();
  
    G4double x = FindScaledEnergy(iz1, rand, logTkin, yymin, yymax);
    if(iz1 != iz2) {
      G4double x2 = FindScaledEnergy(iz2, rand, logTkin, yymin, yymax);
      G4double lz1= nist->GetLOGZ(iz1);
      G4double lz2= nist->GetLOGZ(iz2);
      //G4cout << count << ".  x= " << x << "  x2= " << x2 
      //             << " Z1= " << iz1 << " Z2= " << iz2 << G4endl;
      x += (x2 - x)*(lnZ - lz1)/(lz2 - lz1);
    }
    //G4cout << "x= " << x << "  coeff= " << coeff << G4endl;
    pairEnergy = kinEnergy*G4Exp(x*coeff);
    
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while((pairEnergy < minEnergy || pairEnergy > maxEnergy) && 10 > count);

  //G4cout << "## pairEnergy(GeV)= " << pairEnergy/GeV 
  //         << " Etot(GeV)= " << totalEnergy/GeV << G4endl; 

  // sample r=(mu+mu-)/pairEnergy  ( uniformly .....)
  G4double rmax = 1 - 2*muonMass/(pairEnergy);
  G4double r = rmax * (-1.+2.*G4UniformRand()) ;

  // compute energies from pairEnergy,r
  G4double muMinusEnergy = (1.-r)*pairEnergy*0.5;
  G4double muPlusEnergy = pairEnergy - muMinusEnergy;

  // Sample angles 
  G4ThreeVector muMinusDirection, muPlusDirection;
  //
  GetAngularDistribution()->SamplePairDirections(aDynamicParticle, 
                                                 muMinusEnergy, muPlusEnergy,
                                                 muMinusDirection, muPlusDirection);
  // create G4DynamicParticle object for mu+mu-
  muMinusEnergy = std::max(muMinusEnergy - muonMass, 0.0);
  muPlusEnergy = std::max(muPlusEnergy - muonMass, 0.0);
  G4DynamicParticle* aParticle1 =
    new G4DynamicParticle(theMuonMinus,muMinusDirection,muMinusEnergy);
  G4DynamicParticle* aParticle2 = 
    new G4DynamicParticle(theMuonPlus,muPlusDirection,muPlusEnergy);
  // Fill output vector
  vdp->push_back(aParticle1);
  vdp->push_back(aParticle2);

  // primary change
  kinEnergy -= pairEnergy;
  partDirection *= totalMomentum;
  partDirection -= (aParticle1->GetMomentum() + aParticle2->GetMomentum());
  partDirection = partDirection.unit();

  // if energy transfer is higher than threshold (very high by default)
  // then stop tracking the primary particle and create a new secondary
  if (pairEnergy > SecondaryThreshold()) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    G4DynamicParticle* newdp = 
      new G4DynamicParticle(particle, partDirection, kinEnergy);
    vdp->push_back(newdp);
  } else { // continue tracking the primary e-/e+ otherwise
    fParticleChange->SetProposedMomentumDirection(partDirection);
    fParticleChange->SetProposedKineticEnergy(kinEnergy);
  }
  //G4cout << "--- G4MuonToMuonPairProductionModel::SampleSecondaries done" << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuonToMuonPairProductionModel::DataCorrupted(G4int Z, G4double logTkin) const
{
  G4ExceptionDescription ed;
  ed << "G4ElementData is not properly initialized Z= " << Z
     << " Ekin(MeV)= " << G4Exp(logTkin)
     << " IsMasterThread= " << IsMaster() 
     << " Model " << GetName();
  G4Exception("G4MuonToMuonPairProductionModel","em0033",FatalException,ed,"");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
