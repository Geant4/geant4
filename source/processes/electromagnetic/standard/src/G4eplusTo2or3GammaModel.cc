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
// File name:   G4eplusTo2or3GammaModel
//
// Author:      Vladimir Ivanchenko and Omrame Kadri
//
// Creation date: 29.03.2018
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4eplusTo2or3GammaModel.hh"
#include "G4eplusTo3GammaOKVIModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4TrackStatus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4eplusTo2or3GammaModel::fCrossSection   = nullptr;
G4PhysicsVector* G4eplusTo2or3GammaModel::f3GProbability  = nullptr;

G4eplusTo2or3GammaModel::G4eplusTo2or3GammaModel()
  : G4VEmModel("eplusTo2or3gamma"),
    fDeltaMin(0.001),
    fDelta(fDeltaMin),
    fGammaTh(CLHEP::MeV)
{
  theGamma = G4Gamma::Gamma();
  fParticleChange = nullptr;
  f3GModel = new G4eplusTo3GammaOKVIModel();
  SetTripletModel(f3GModel);

  // instantiate vectors once
  if (nullptr == fCrossSection) {
    G4double emin = 10*CLHEP::eV;
    G4double emax = 100*CLHEP::TeV;
    G4int nbins = 20*G4lrint(std::log10(emax/emin));
    fCrossSection = new G4PhysicsLogVector(emin, emax, nbins, true);
    f3GProbability = new G4PhysicsLogVector(emin, emax, nbins, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusTo2or3GammaModel::~G4eplusTo2or3GammaModel()
{
  if (IsMaster()) {
    delete fCrossSection;
    delete f3GProbability;
    fCrossSection = nullptr;
    f3GProbability = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusTo2or3GammaModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& cuts)
{
  // here particle change is set for the triplet model
  if (nullptr == fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }
  // initialialise 3-gamma model before new run
  f3GModel->Initialise(p, cuts);
  fGammaTh = G4EmParameters::Instance()->LowestTripletEnergy();

  // initialise vectors
  if (IsMaster()) {
    std::size_t num = fCrossSection->GetVectorLength();
    for (std::size_t i=0; i<num; ++i) {
      G4double e = fCrossSection->Energy(i);
      G4double cs2 = ComputeCrossSectionPerElectron(e);
      G4double cs3 = f3GModel->ComputeCrossSectionPerElectron(e);
      cs2 += cs3;
      fCrossSection->PutValue(i, cs2);
      G4double y = (cs2 > 0.0) ? cs3/cs2 : 0.0;
      f3GProbability->PutValue(i, y);
    }
    fCrossSection->FillSecondDerivatives();
    f3GProbability->FillSecondDerivatives();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4eplusTo2or3GammaModel::ComputeCrossSectionPerElectron(G4double kinEnergy)
{
  // Calculates the cross section per electron of annihilation into two 
  // photons from the Heilter formula with the radiation correction to 3 gamma 
  // annihilation channel. (A.A.) rho is changed

  G4double ekin   = std::max(CLHEP::eV, kinEnergy);   
  G4double tau    = ekin/CLHEP::electron_mass_c2;
  G4double gam    = tau + 1.0;
  G4double gamma2 = gam*gam;
  G4double bg2    = tau * (tau+2.0);
  G4double bg     = std::sqrt(bg2);
  G4double rho = (gamma2+4.*gam+1.)*G4Log(gam+bg)/(gamma2-1.) 
    - (gam+3.)/(std::sqrt(gam*gam - 1.));
  G4double eGammaCMS = CLHEP::electron_mass_c2 * std::sqrt(0.5*(tau + 2.0));
  fDelta = std::max(fDeltaMin, fGammaTh/eGammaCMS);
  f3GModel->SetDelta(fDelta);

  static const G4double pir2 =
    CLHEP::pi*CLHEP::classic_electr_radius*CLHEP::classic_electr_radius;
  G4double cross = (pir2*rho + alpha_rcl2*2.*G4Log(fDelta)*rho*rho)/(gam+1.);

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eplusTo2or3GammaModel::ComputeCrossSectionPerAtom(
                                    const G4ParticleDefinition*,
                                    G4double kineticEnergy, G4double Z,
				    G4double, G4double, G4double)
{
  // Calculates the cross section per atom of annihilation into two photons
  G4double cross = Z*fCrossSection->Value(kineticEnergy);
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eplusTo2or3GammaModel::CrossSectionPerVolume(
					const G4Material* material,
					const G4ParticleDefinition*,
					      G4double kineticEnergy,
					      G4double, G4double)
{
  // Calculates the cross section per volume of annihilation into two photons
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*fCrossSection->Value(kineticEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward, 
// Nature 4065 (1947) 435.

void G4eplusTo2or3GammaModel::SampleSecondaries(
				     std::vector<G4DynamicParticle*>* vdp,
                                     const G4MaterialCutsCouple* couple,
                                     const G4DynamicParticle* dp,
                                     G4double, G4double)
{
  // kill primary positron
  fParticleChange->SetProposedKineticEnergy(0.0);
  fParticleChange->ProposeTrackStatus(fStopAndKill);

  G4ThreeVector posiDirection = dp->GetMomentumDirection();

  // Case at rest not considered anymore
  G4double posiKinEnergy = dp->GetKineticEnergy();
  G4LorentzVector lv(dp->GetMomentum(),
                     posiKinEnergy + 2*CLHEP::electron_mass_c2);
  G4double eGammaCMS = 0.5 * lv.mag();

  // 3-gamma annihilation
  if (G4UniformRand() < f3GProbability->Value(posiKinEnergy)) {
    fDelta = std::max(fDeltaMin, fGammaTh/eGammaCMS);
    f3GModel->SetDelta(fDelta);
    f3GModel->SampleSecondaries(vdp, couple, dp);
    return;
  }

  // 2-gamma annihilation
  G4double tau     = posiKinEnergy/CLHEP::electron_mass_c2;
  G4double gam     = tau + 1.0;
  G4double tau2    = tau + 2.0;
  G4double sqgrate = std::sqrt(tau/tau2)*0.5;
  G4double sqg2m1  = std::sqrt(tau*tau2);

  // limits of the energy sampling
  G4double epsilmin = 0.5 - sqgrate;
  G4double epsilmax = 0.5 + sqgrate;
  G4double epsilqot = epsilmax/epsilmin;

  //
  // sample the energy rate of the created gammas
  //
  G4double epsil, greject;

  do {
    epsil = epsilmin*G4Exp(G4Log(epsilqot)*G4UniformRand());
    greject = 1. - epsil + (2.*gam*epsil-1.)/(epsil*tau2*tau2);
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( greject < G4UniformRand());

  //
  // scattered Gamma angles. ( Z - axis along the parent positron)
  //
  G4double cost = (epsil*tau2-1.)/(epsil*sqg2m1);
  if (std::abs(cost) > 1.0) {
    G4cout << "### G4eeToTwoGammaModel WARNING cost= " << cost
	   << " positron Ekin(MeV)= " << posiKinEnergy
	   << " gamma epsil= " << epsil
	   << G4endl;
    if (cost > 1.0) cost = 1.0;
    else cost = -1.0; 
  }
  G4double sint = std::sqrt((1.+cost)*(1.-cost));
  G4double phi  = CLHEP::twopi * G4UniformRand();

  //
  // kinematic of the created pair
  //
  G4double totalEnergy = posiKinEnergy + 2.0*CLHEP::electron_mass_c2;
  G4double phot1Energy = epsil*totalEnergy;

  G4ThreeVector phot1Direction(sint*std::cos(phi), sint*std::sin(phi), cost);
  phot1Direction.rotateUz(posiDirection);
  auto aGamma1 = new G4DynamicParticle (theGamma,phot1Direction, phot1Energy);
  phi = CLHEP::twopi * G4UniformRand();
  G4double cosphi = std::cos(phi);
  G4double sinphi = std::sin(phi);
  G4ThreeVector pol(cosphi, sinphi, 0.0);
  pol.rotateUz(phot1Direction);
  aGamma1->SetPolarization(pol.x(),pol.y(),pol.z());

  G4double phot2Energy = (1.-epsil)*totalEnergy;
  G4double posiP = std::sqrt(posiKinEnergy*(posiKinEnergy+2.*electron_mass_c2));
  G4ThreeVector dir = posiDirection*posiP - phot1Direction*phot1Energy;
  G4ThreeVector phot2Direction = dir.unit();

  // create G4DynamicParticle object for the particle2
  auto aGamma2 = new G4DynamicParticle (theGamma, phot2Direction, phot2Energy);

  //!!! likely problematic direction to be checked
  pol.set(-sinphi, cosphi, 0.0);
  pol.rotateUz(phot1Direction);
  cost = pol*phot2Direction;
  pol -= cost*phot2Direction;
  pol = pol.unit();
  aGamma2->SetPolarization(pol.x(),pol.y(),pol.z());
  /*                                                                                                                                                                         
    G4cout << "Annihilation on fly: e0= " << posiKinEnergy
           << " m= " << electron_mass_c2
           << " e1= " << phot1Energy
           << " e2= " << phot2Energy << " dir= " << dir                                                                                                                      
           << " -> " << phot1Direction << " "
           << phot2Direction << G4endl;
    */

  vdp->push_back(aGamma1);
  vdp->push_back(aGamma2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
