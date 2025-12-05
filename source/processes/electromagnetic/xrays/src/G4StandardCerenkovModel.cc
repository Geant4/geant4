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
// G4StandardCerenkovModel
//
// Created 25.05.2025 V.Ivanchenko
//
// --------------------------------------------------------------------

#include "G4StandardCerenkovModel.hh"

#include "G4ios.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4Poisson.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4OpticalPhoton.hh"

namespace
{
  constexpr G4double minAllowedStep = 0.001*CLHEP::mm;
  constexpr G4int nvec = 10; // number of slices in beta of projectile
}

std::vector<G4double>* G4StandardCerenkovModel::fBetaLim = nullptr;
std::vector<G4MaterialPropertyVector*>* G4StandardCerenkovModel::fRindex = nullptr;
std::vector<std::vector<G4double>* >* G4StandardCerenkovModel::fMeanNumberOfPhotons = nullptr;
std::vector<std::vector<std::vector<G4double>* >* >* G4StandardCerenkovModel::fIntegral = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StandardCerenkovModel::G4StandardCerenkovModel()
  : G4VXRayModel("theCerenkov")
{
  if (nullptr == fBetaLim) {
    isInitializer = true;
    fBetaLim = new std::vector<G4double>;
    fMeanNumberOfPhotons = new std::vector<std::vector<G4double>* >;
    fIntegral = new std::vector<std::vector<std::vector<G4double>* >* >;
    fRindex = new std::vector<G4MaterialPropertyVector*>;
  }
  fPhoton = G4OpticalPhoton::OpticalPhoton();
  fRfact = 369.81 / (CLHEP::eV * CLHEP::cm); // number of photons per mm
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StandardCerenkovModel::~G4StandardCerenkovModel()
{
  if (isInitializer && isInitialized) {
    delete fBetaLim;
    for (auto const & ptr : *fRindex) {
      delete ptr;
    }
    delete fRindex;
    for (auto const & ptr : *fMeanNumberOfPhotons) {
      delete ptr;
    }
    delete fMeanNumberOfPhotons;
    for (auto const & ptr : *fIntegral) {
      for (auto const & p : *ptr) {
	delete p;
      }
      delete ptr;
    }
    delete fIntegral;
    fIntegral = nullptr;
    fBetaLim = nullptr;
    fRindex = nullptr;
    fMeanNumberOfPhotons = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4StandardCerenkovModel::InitialiseModel()
{
  G4double beta = 1.0;
  if (0 == nVolumes) { return; }
  if (isInitializer && !isInitialized) {
    isInitialized = true;
    fBetaLim->resize(nVolumes, 1.0);
    fRindex->resize(nVolumes, nullptr);
    fMeanNumberOfPhotons->resize(nVolumes, new std::vector<G4double>(nvec, 0.0));
    fIntegral->resize(nVolumes, nullptr);
      
    for (std::size_t i = 0; i < nVolumes; ++i) {
      auto mat = ((*pLogicalVolumes)[i])->GetMaterial();
      auto MPT = mat->GetMaterialPropertiesTable();
	
      if (nullptr == MPT) { continue; }
      G4MaterialPropertyVector* rindex = MPT->GetProperty(kRINDEX);
      if (nullptr == rindex) { continue; }
      (*fRindex)[i] = rindex;

      G4double nMax = rindex->GetMaxValue();
      if (nMax <= 1.0) { continue; }
      G4double b = 1.0/nMax;
      (*fBetaLim)[i] = b;
      beta = std::min(beta, b);
      G4double dbeta = (1.0 - b)/(G4double)(nvec - 1);
      G4double b0 = b;
      std::size_t nn = rindex->GetVectorLength();
      if (0 == nn) { continue; }

      auto iptr = new std::vector<std::vector<G4double>* >((std::size_t)nvec, nullptr);
      (*fIntegral)[i] = iptr;
      for (auto & ptr : *iptr) {
	ptr =  new std::vector<G4double>(nn, 0.0);
      }
      
      // Initialisation for charge = 1.0
      G4double y0 = AverageNumberOfPhotons(1.0, b0, (*rindex)[0]);
      (*(*fMeanNumberOfPhotons)[i])[0] = y0;
      if (1 == nn) { continue; }
      
      // Initialisation for charge = 1.0
      for (G4int k = 0; k < nvec; ++k) {
	(*((*fIntegral)[i]))[k]->resize(nn, 0.0);
	G4double sum = 0.0;
	G4double e0 = rindex->GetMinEnergy();
	G4double deltae = rindex->GetMaxEnergy() - e0;
	for (std::size_t j = 1; j < nn; ++j) {
	  G4double e = rindex->Energy(j);
	  G4double y = AverageNumberOfPhotons(1.0, beta, (*rindex)[j]);
	  sum += 0.5*(y - y0)*(e - e0);
	  y0 = y;
	  e0 = e;
	  (*(*((*fIntegral)[i]))[k])[j] = sum;
	}
	if (deltae > 0.0) {  (*(*fMeanNumberOfPhotons)[i])[k] = sum/deltae; }
	if (sum > 0.0) { sum = 1.0/sum; }
	for (std::size_t j = 1; j < nn; ++j) {
	  G4double y = (*(*((*fIntegral)[i]))[k])[j]/sum;
	  (*(*((*fIntegral)[i]))[k])[j] = y;
	}
	beta += dbeta;
	beta = std::min(beta, 1.0);
      }
    }
  } else {
    for (std::size_t i = 0; i < nVolumes; ++i) {
      beta = std::min(beta, (*fBetaLim)[i]);
    }
  }
  pBetaMin = beta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4StandardCerenkovModel::StepLimitForVolume(G4double& limit)
{
  G4double betaMin = (*fBetaLim)[pIndex];
  if (pPreStepBeta <= betaMin) { return false; }

  // step limitation 
  betaMin = std::max(betaMin, pPreStepBeta*pMaxBetaChange);
  
  auto const dynPart = pCurrentTrack->GetDynamicParticle();
  fParticle = dynPart->GetDefinition();
  fPreStepKinE = dynPart->GetKineticEnergy();
  fCharge = fParticle->GetPDGCharge()/CLHEP::eplus;
  fMass = fParticle->GetPDGMass();
  G4double x = limit;

  // If the step is smaller than G4ThreeVector::getTolerance(), it may happen
  // that the particle does not move. See bug 1992.
  if (x < minAllowedStep) { return false; }

  // If user has defined an average maximum number of photons to be generated in
  // a Step, then calculate the Step length for that number of photons using preStep beta.
  G4double nmax = (*fRindex)[pIndex]->GetMaxValue();
  fMeanNPhotons = x * AverageNumberOfPhotons(fCharge, pPreStepBeta, nmax);
  if (pMaxPhotons < fMeanNPhotons) {
    x *= pMaxPhotons/fMeanNPhotons;
    x = std::max(x, minAllowedStep);
  }
  limit = x;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4StandardCerenkovModel::SampleXRays(std::vector<G4Track*>& out,
					  const G4Step& step)
{
  auto const preStep = step.GetPreStepPoint();
  auto const postStep = step.GetPostStepPoint();
  G4double kinE = postStep->GetKineticEnergy();
  G4double beta = pPreStepBeta;
  G4double b = (*fBetaLim)[pIndex];
  G4double dbeta = (1.0 - b)/(G4double)(nvec - 1);
  G4int idx = std::max(G4int((beta - b)/dbeta), 0);
  idx = std::min(idx, nvec - 1); 
  std::vector<G4double>* v = (*((*fIntegral)[pIndex]))[idx];
  
  G4ThreeVector pos = preStep->GetPosition();
  G4ThreeVector dpos = postStep->GetPosition() - pos;
  G4ThreeVector dir = dpos.unit();
  G4double time = preStep->GetGlobalTime();
  G4double dt = postStep->GetGlobalTime() - time;
  G4double de = fPreStepKinE - kinE;

  // define sub-steps inside the current step
  G4double x = (1.0 + fPreStepKinE/fMass);
  G4double delim = fMass*pMaxBetaChange*x*x*x;
  G4int nn = (G4int)(de/delim) + 1;
  x = 1.0/(G4double)nn;
  dpos *= x;
  dt *= x;
  de *= x;
  G4double delta = dpos.mag();
  fMeanNPhotons *= x;
  
  // produce Cerenkov gamma - loop over sub-steps
  G4double ekin = fPreStepKinE;
  auto const rindex = (*fRindex)[pIndex];
  std::size_t ni = rindex->GetVectorLength();
  if (0 == ni) { return; }
  G4double emin = rindex->Energy(0);
  G4double mean = delta*fCharge*fCharge*(*(*fMeanNumberOfPhotons)[pIndex])[idx];
  
  for (G4int i=0; i<nn; ++i) {
    G4int ngamma = (G4int)G4Poisson(mean);
    for (G4int j = 0; j < ngamma; ++j) {
      G4double q = G4UniformRand();
      G4double e = emin;
      G4double n = (*rindex)[0];
      if (ni > 1) {
	for (std::size_t k = 1; k < ni; ++k) {
	  if ((*v)[k] <= q) {
	    e = rindex->Energy(k - 1);
	    e += (rindex->Energy(k) - e)*(q - (*v)[k - 1])/((*v)[k] - (*v)[k - 1]);
	    n = (*rindex)[k - 1];
	    n += ((*rindex)[k] - n)*(q - (*v)[k - 1])/((*v)[k] - (*v)[k - 1]);
	  }
	}
      }
      q = G4UniformRand();
      G4double t = time + q*dt;
      G4ThreeVector posnew = pos + q*dpos;
      G4double minCos = 1.0/(n*beta);
      G4double maxSin2 = (1.0 - minCos)*(1.0 + minCos);
      G4double cost, sint2;
      do {
        cost = 1.0 - G4UniformRand()*(1.0 - minCos);
	sint2 = (1.0 - cost)*(1.0 + cost);
	q = G4UniformRand();
      } while(q * maxSin2 > sint2);
      
      G4double sint = std::sqrt(sint2);
      G4double phi = G4UniformRand()*CLHEP::twopi;
      G4double cosPhi = std::cos(phi);
      G4double sinPhi = std::sin(phi);
      G4ThreeVector dirnew(sint*cosPhi, sint*sinPhi, cost);
      dirnew.rotateUz(dir);

      // Determine polarization of new photon
      G4ThreeVector photonPolarization(cost*cosPhi, cost*sinPhi, -sint);

      // Rotate back to original coord system
      photonPolarization.rotateUz(dir);

      auto dp = new G4DynamicParticle(fPhoton, dirnew, e);
      dp->SetPolarization(photonPolarization);
      auto track = new G4Track(dp, t, posnew);
      out.push_back(track);
    }
    pos += dpos;
    time += dt;
    ekin -= de;
    beta = std::sqrt(ekin * (ekin + 2*fMass))/(fMass + ekin);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4StandardCerenkovModel::ModelDescription(std::ostream& out) const
{
  out << "The Cerenkov effect simulates optical photons created by the\n";
  out << "passage of charged particles through matter. Materials need\n";
  out << "to have the property RINDEX (refractive index) defined." << G4endl;
}




