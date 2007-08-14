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
// $Id: G4mplIonisationModel.cc,v 1.4 2007-08-14 09:47:58 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4mplIonisationModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 06.09.2005
//
// Modifications:
// 12.08.2007 Changing low energy approximation and extrapolation. Small bug fixing and refactoring (M. Vladymyrov)
//
//
// -------------------------------------------------------------------
// References
// [1] Steven P. Ahlen: Energy loss of relativistic heavy ionizing particles, Rev. Mod. Phys 52(1980), p121



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4mplIonisationModel.hh"
#include "Randomize.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForLoss.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4mplIonisationModel::G4mplIonisationModel(G4double mCharge, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
  magCharge(mCharge),
  twoln10(log(100.0)),
  beta2low(0.004),
  beta2lim(0.02),
  bg2lim(beta2lim*(1.0 + beta2lim)),
  lastMatelial(0)
{
  nmpl         = G4int(abs(magCharge) * 2 * fine_structure_const + 0.5);
  if(nmpl > 6)      nmpl = 6;
  else if(nmpl < 1) nmpl = 1;
  pi_hbarc2_over_mc2 = pi * hbarc * hbarc / electron_mass_c2;
  chargeSquare = magCharge * magCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4mplIonisationModel::~G4mplIonisationModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisationModel::Initialise(const G4ParticleDefinition* p,
				      const G4DataVector&)
{
  monopole = p;
  mass     = monopole->GetPDGMass();

  if(pParticleChange) 
    fParticleChange = reinterpret_cast<G4ParticleChangeForLoss*>(pParticleChange);
  else 
    fParticleChange = new G4ParticleChangeForLoss();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::ComputeDEDXPerVolume(const G4Material* material,
						    const G4ParticleDefinition*,
						    G4double kineticEnergy,
						    G4double)
{
  G4double tau   = kineticEnergy / mass;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau + 2.0);
  G4double beta2 = bg2 / (gam * gam);

  // approximation in low energy region
  G4double dedx0 = 0, dedx = 0;

  if(beta2 > beta2low) {
    dedx = ComputeDEDXAhlen(material, bg2);
  }

  if(beta2 < beta2lim) {
    if(lastMatelial != material){
      approxConst = ComputeDEDXAhlen(material, bg2lim);
      lastMatelial = (G4Material *)material;
    };
    dedx = dedx0 = approxConst * sqrt(sqrt(abs(beta2 / bg2lim)));

    // extrapolation between two formula 
    if(beta2 > beta2low) {
   
      G4double kapa1 = sqrt(abs(beta2)), kapa2 = kapa1;
      kapa1 -= sqrt(beta2low);
      kapa2 -= sqrt(beta2lim);
      kapa1 *= kapa1;
      kapa2 *= kapa2;

      G4double x = kapa1 * dedx + kapa2 * dedx0;
      x /= kapa1 + kapa2;
      dedx = x;
    }
  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::ComputeDEDXAhlen(const G4Material* material, G4double bg2)
{

  G4double eDensity = material->GetElectronDensity();
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double cden  = material->GetIonisation()->GetCdensity();
  G4double mden  = material->GetIonisation()->GetMdensity();
  G4double aden  = material->GetIonisation()->GetAdensity();
  G4double x0den = material->GetIonisation()->GetX0density();
  G4double x1den = material->GetIonisation()->GetX1density();

  // Ahlen's formula for nonconductors, [1]p157, f(5.7)
  G4double dedx = log(2.0 * electron_mass_c2 * bg2 / eexc) - 0.5;

  // Kazama et al. cross-section correction
  G4double  k = 0.406;
  if(nmpl > 1) k = 0.346;

  // Bloch correction
  const G4double B[7] = { 0.0, 0.248, 0.672, 1.022, 1.243, 1.464, 1.685}; 

  dedx += 0.5 * k - B[nmpl];

  // density effect correction
  G4double deltam;
  G4double x = log(bg2) / twoln10;
  if ( x >= x0den ) {
    deltam = twoln10 * x - cden;
    if ( x < x1den ) deltam += aden * pow((x1den-x), mden);
    dedx -= 0.5 * deltam;
  }

  // now compute the total ionization loss
  dedx *=  pi_hbarc2_over_mc2 * eDensity * nmpl * nmpl;

  if (dedx < 0.0) dedx = 0;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::SampleFluctuations(
				       const G4Material* material,
				       const G4DynamicParticle* dp,
				       G4double& tmax,
				       G4double& length,
				       G4double& meanLoss)
{
  G4double siga = Dispersion(material,dp,tmax,length);
  G4double loss = meanLoss;
  siga = sqrt(siga);
  G4double twomeanLoss = meanLoss + meanLoss;

  if(twomeanLoss < siga) {
    G4double x;
    do {
      loss = twomeanLoss*G4UniformRand();
      x = (loss - meanLoss)/siga;
    } while (1.0 - 0.5*x*x < G4UniformRand());
  } else {
    do {
      loss = G4RandGauss::shoot(meanLoss,siga);
    } while (0.0 > loss || loss > twomeanLoss);
  }
  return loss;
}
