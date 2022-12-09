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
// File name:     G4mplIonisationModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 06.09.2005
//
// Modifications:
// 12.08.2007 Changing low energy approximation and extrapolation. 
//            Small bug fixing and refactoring (M. Vladymyrov)
// 13.11.2007 Use low-energy asymptotic from [3] (V.Ivanchenko) 
//
//
// -------------------------------------------------------------------
// References
// [1] Steven P. Ahlen: Energy loss of relativistic heavy ionizing particles, 
//     S.P. Ahlen, Rev. Mod. Phys 52(1980), p121
// [2] K.A. Milton arXiv:hep-ex/0602040
// [3] S.P. Ahlen and K. Kinoshita, Phys. Rev. D26 (1982) 2347


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4mplIonisationModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4double>* G4mplIonisationModel::dedx0 = nullptr;

G4mplIonisationModel::G4mplIonisationModel(G4double mCharge, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
  magCharge(mCharge),
  twoln10(G4Log(100.0)),
  betalow(0.01),
  betalim(0.1),
  beta2lim(betalim*betalim),
  bg2lim(beta2lim*(1.0 + beta2lim))
{
  nmpl = G4int(std::abs(magCharge) * 2 * CLHEP::fine_structure_const + 0.5);
  if(nmpl > 6)      { nmpl = 6; }
  else if(nmpl < 1) { nmpl = 1; }
  pi_hbarc2_over_mc2 = CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc/CLHEP::electron_mass_c2;
  chargeSquare = magCharge * magCharge;
  dedxlim = 45.*nmpl*nmpl*CLHEP::GeV*CLHEP::cm2/CLHEP::g;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4mplIonisationModel::~G4mplIonisationModel()
{
  if(IsMaster()) { delete dedx0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisationModel::SetParticle(const G4ParticleDefinition* p)
{
  monopole = p;
  mass     = monopole->GetPDGMass();
  G4double emin = 
    std::min(LowEnergyLimit(),0.1*mass*(1./std::sqrt(1. - betalow*betalow) - 1.)); 
  G4double emax = 
    std::max(HighEnergyLimit(),10.*mass*(1./std::sqrt(1. - beta2lim) - 1.)); 
  SetLowEnergyLimit(emin);
  SetHighEnergyLimit(emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisationModel::Initialise(const G4ParticleDefinition* p,
				      const G4DataVector&)
{
  if(nullptr == monopole) { SetParticle(p); }
  if(nullptr == fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }
  if(IsMaster()) {
    if(nullptr == dedx0) { dedx0 = new std::vector<G4double>; }
    G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
    G4int n = (G4int)dedx0->size();
    if(n < numOfCouples) { dedx0->resize(numOfCouples); }

    G4Pow* g4calc = G4Pow::GetInstance();

    // initialise vector assuming low conductivity
    for(G4int i=0; i<numOfCouples; ++i) {

      const G4Material* material = 
        theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      G4double eDensity = material->GetElectronDensity();
      G4double vF2 = 2*electron_Compton_length*g4calc->A13(3.*pi*pi*eDensity);
      (*dedx0)[i] = pi_hbarc2_over_mc2*eDensity*nmpl*nmpl*
        (G4Log(vF2/fine_structure_const) - 0.5)/vF2;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::ComputeDEDXPerVolume(const G4Material* material,
						    const G4ParticleDefinition* p,
						    G4double kineticEnergy,
						    G4double)
{
  if(nullptr == monopole) { SetParticle(p); }
  G4double tau   = kineticEnergy / mass;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau + 2.0);
  G4double beta2 = bg2 / (gam * gam);
  G4double beta  = std::sqrt(beta2);

  // low-energy asymptotic formula
  //G4double dedx  = dedxlim*beta*material->GetDensity();
  G4double dedx = (*dedx0)[CurrentCouple()->GetIndex()]*beta;

  // above asymptotic
  if(beta > betalow) {

    // high energy
    if(beta >= betalim) {
      dedx = ComputeDEDXAhlen(material, bg2);

    } else {

      //G4double dedx1 = dedxlim*betalow*material->GetDensity();
      G4double dedx1 = (*dedx0)[CurrentCouple()->GetIndex()]*betalow;
      G4double dedx2 = ComputeDEDXAhlen(material, bg2lim);

      // extrapolation between two formula 
      G4double kapa2 = beta - betalow;
      G4double kapa1 = betalim - beta;
      dedx = (kapa1*dedx1 + kapa2*dedx2)/(kapa1 + kapa2);
    }
  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::ComputeDEDXAhlen(const G4Material* material, 
						G4double bg2)
{
  G4double eDensity = material->GetElectronDensity();
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double cden  = material->GetIonisation()->GetCdensity();
  G4double mden  = material->GetIonisation()->GetMdensity();
  G4double aden  = material->GetIonisation()->GetAdensity();
  G4double x0den = material->GetIonisation()->GetX0density();
  G4double x1den = material->GetIonisation()->GetX1density();

  // Ahlen's formula for nonconductors, [1]p157, f(5.7)
  G4double dedx = std::log(2.0 * electron_mass_c2 * bg2 / eexc) - 0.5;

  // Kazama et al. cross-section correction
  G4double  k = 0.406;
  if(nmpl > 1) k = 0.346;

  // Bloch correction
  const G4double B[7] = { 0.0, 0.248, 0.672, 1.022, 1.243, 1.464, 1.685}; 

  dedx += 0.5 * k - B[nmpl];

  // density effect correction
  G4double deltam;
  G4double x = std::log(bg2) / twoln10;
  if ( x >= x0den ) {
    deltam = twoln10 * x - cden;
    if ( x < x1den ) deltam += aden * std::pow((x1den-x), mden);
    dedx -= 0.5 * deltam;
  }

  // now compute the total ionization loss
  dedx *=  pi_hbarc2_over_mc2 * eDensity * nmpl * nmpl;

  if (dedx < 0.0) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
					     const G4MaterialCutsCouple*,
					     const G4DynamicParticle*,
					     G4double,
					     G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::SampleFluctuations(
				       const G4MaterialCutsCouple* couple,
				       const G4DynamicParticle* dp,
                                       const G4double tcut,
                                       const G4double tmax,
				       const G4double length,
				       const G4double meanLoss)
{
  G4double siga = Dispersion(couple->GetMaterial(),dp,tcut,tmax,length);
  G4double loss = meanLoss;
  siga = std::sqrt(siga);
  G4double twomeanLoss = meanLoss + meanLoss;

  if(twomeanLoss < siga) {
    G4double x;
    do {
      loss = twomeanLoss*G4UniformRand();
      x = (loss - meanLoss)/siga;
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (1.0 - 0.5*x*x < G4UniformRand());
  } else {
    do {
      loss = G4RandGauss::shoot(meanLoss,siga);
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (0.0 > loss || loss > twomeanLoss);
  }
  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4mplIonisationModel::Dispersion(const G4Material* material,
					  const G4DynamicParticle* dp,
					  const G4double tcut,
					  const G4double tmax,
					  const G4double length)
{
  G4double siga = 0.0;
  G4double tau   = dp->GetKineticEnergy()/mass;
  if(tau > 0.0) { 
    const G4double beta = dp->GetBeta();
    siga  = (tmax/(beta*beta) - 0.5*tcut) * twopi_mc2_rcl2 * length
      * material->GetElectronDensity() * chargeSquare;
  }
  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
