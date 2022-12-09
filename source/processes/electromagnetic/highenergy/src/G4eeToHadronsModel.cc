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
// File name:     G4eeToHadronsModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 12.08.2003
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-05-05 Use optimized interfaces (V.Ivantchenko)
//
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToHadronsModel.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4PionPlus.hh"
#include "Randomize.hh"
#include "G4Vee2hadrons.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToHadronsModel::G4eeToHadronsModel(G4Vee2hadrons* mod, G4int ver,
                                       const G4String& nam)
  : G4VEmModel(nam),
    model(mod),
    verbose(ver)
{
  theGamma = G4Gamma::Gamma();
  highKinEnergy = HighEnergyLimit();
  lowKinEnergy  = LowEnergyLimit();
  emin = lowKinEnergy;
  emax = highKinEnergy;
  peakKinEnergy = highKinEnergy;
  epeak = emax;
  //verbose = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadronsModel::~G4eeToHadronsModel()
{
  delete model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::Initialise(const G4ParticleDefinition*,
                                    const G4DataVector&)
{
  if(isInitialised) { return; }
  isInitialised  = true;

  // CM system
  emin = model->LowEnergy();
  emax = model->HighEnergy();

  // peak energy
  epeak = std::min(model->PeakEnergy(), emax);
  
  if(verbose>0) {
    G4cout << "G4eeToHadronsModel::Initialise: " << G4endl;
    G4cout << "CM System: emin(MeV)= " << emin/MeV
           << " epeak(MeV)= " << epeak/MeV
           << " emax(MeV)= " << emax/MeV
           << G4endl;
  }

  crossBornPerElectron = model->PhysicsVector();
  crossPerElectron     = model->PhysicsVector(); 
  nbins = (G4int)crossPerElectron->GetVectorLength();
  for(G4int i=0; i<nbins; ++i) {
    G4double e  = crossPerElectron->Energy(i);
    G4double cs = model->ComputeCrossSection(e);
    crossBornPerElectron->PutValue(i, cs);
  }
  ComputeCMCrossSectionPerElectron();

  if(verbose>1) {
    G4cout << "G4eeToHadronsModel: Cross sections per electron"
           << " nbins= " << nbins
           << " emin(MeV)= " << emin/MeV
           << " emax(MeV)= " << emax/MeV
           << G4endl;
    for(G4int i=0; i<nbins; ++i) {
      G4double e  = crossPerElectron->Energy(i);
      G4double s1 = crossPerElectron->Value(e);
      G4double s2 = crossBornPerElectron->Value(e);
      G4cout << "E(MeV)= " << e/MeV
             << "  cross(nb)= " << s1/nanobarn
             << "  crossBorn(nb)= " << s2/nanobarn
	     << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::CrossSectionPerVolume(
				      const G4Material* mat,
				      const G4ParticleDefinition* p,
				      G4double kineticEnergy,
				      G4double, G4double)
{
  return mat->GetElectronDensity()*
    ComputeCrossSectionPerElectron(p, kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::ComputeCrossSectionPerAtom(
                                      const G4ParticleDefinition* p,
				      G4double kineticEnergy,
				      G4double Z, G4double,
				      G4double, G4double)
{
  return Z*ComputeCrossSectionPerElectron(p, kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::ComputeCrossSectionPerElectron(
                                          const G4ParticleDefinition*,
                                                G4double energy,
                                                G4double, G4double)
{
  return crossPerElectron->Value(energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
					   const G4MaterialCutsCouple*,
					   const G4DynamicParticle* dParticle,
					   G4double,
					   G4double)
{
  G4double t = dParticle->GetKineticEnergy() + 2*electron_mass_c2;
  G4LorentzVector inlv = dParticle->Get4Momentum() + 
    G4LorentzVector(0.0,0.0,0.0,electron_mass_c2);
  G4double e = inlv.m();
  G4ThreeVector inBoost = inlv.boostVector();
  //G4cout << "G4eeToHadronsModel::SampleSecondaries e= " << e 
  //	   << " " << inlv << " " << inBoost <<G4endl;
  if(e > emin) {
    G4DynamicParticle* gamma = GenerateCMPhoton(e);
    G4LorentzVector gLv = gamma->Get4Momentum();
    G4LorentzVector lv(0.0,0.0,0.0,e);
    lv -= gLv;
    G4double mass = lv.m();
    //G4cout << "mass= " << mass << " " << lv << G4endl;
    G4ThreeVector boost = lv.boostVector();
    //G4cout << "mass= " << mass << " " << boost << G4endl;
    const G4ThreeVector dir = gamma->GetMomentumDirection();
    model->SampleSecondaries(newp, mass, dir);
    std::size_t np = newp->size();
    for(std::size_t j=0; j<np; ++j) {
      G4DynamicParticle* dp = (*newp)[j];
      G4LorentzVector v = dp->Get4Momentum();
      v.boost(boost);
      //G4cout << j << ". " << v << G4endl;
      v.boost(inBoost);
      //G4cout << "   " << v << G4endl;
      dp->Set4Momentum(v);
      t -= v.e();
    }
    //G4cout << "Gamma   " << gLv << G4endl;
    gLv.boost(inBoost);
    //G4cout << "        " << gLv << G4endl;
    gamma->Set4Momentum(gLv);
    t -= gLv.e();
    newp->push_back(gamma);
    if(std::abs(t) > CLHEP::MeV) {
      G4cout << "G4eeToHadronsModel::SampleSecondaries: Ebalance(MeV)= " 
	     << t/MeV << " primary 4-momentum: " << inlv <<  G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::ComputeCMCrossSectionPerElectron()
{
  for(G4int i=0; i<nbins; i++) {
    G4double e  = crossPerElectron->Energy(i);
    G4double cs = 0.0;
    if(i > 0) {
      G4double LL   = 2.0*G4Log(e/electron_mass_c2);
      G4double bt  = 2.0*fine_structure_const*(LL - 1.0)/pi;
      G4double btm1= bt - 1.0;
      G4double del = 1. + fine_structure_const*(1.5*LL + pi*pi/3. -2.)/pi;
      G4double s1  = crossBornPerElectron->Value(e);
      G4double e1  = crossPerElectron->Energy(i-1);
      G4double x1  = 1. - e1/e;
      cs += s1*(del*G4Exp(G4Log(x1)*bt) - bt*(x1 - 0.25*x1*x1));
      if(i > 1) {
	G4double e2  = e1;
	G4double x2  = x1;
	G4double s2  = crossBornPerElectron->Value(e2);
	G4double w2  = bt*(del*G4Exp(G4Log(x2)*btm1) - 1.0 + 0.5*x2);
        G4double w1;      

	for(G4int j=i-2; j>=0; --j) {
	  e1  = crossPerElectron->Energy(j);
	  x1  = 1. - e1/e;
	  s1  = crossBornPerElectron->Value(e1);
	  w1  = bt*(del*G4Exp(G4Log(x1)*btm1) - 1.0 + 0.5*x1);
	  cs += 0.5*(x1 - x2)*(w2*s2 + w1*s1);
	  e2 = e1;
	  x2 = x1;
	  s2 = s1;
	  w2 = w1;
	}
      }
    }
    crossPerElectron->PutValue(i, cs);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4eeToHadronsModel::GenerateCMPhoton(G4double e)
{
  G4double x;
  G4DynamicParticle* gamma = nullptr;
  G4double LL   = 2.0*G4Log(e/electron_mass_c2);
  G4double bt  = 2.0*fine_structure_const*(LL - 1.)/pi;
  G4double btm1= bt - 1.0;
  G4double del = 1. + fine_structure_const*(1.5*LL + pi*pi/3. -2.)/pi;

  G4double s0 = crossBornPerElectron->Value(e);
  G4double de = (emax - emin)/(G4double)nbins;
  G4double xmax = 0.5*(1.0 - (emin*emin)/(e*e));
  G4double xmin = std::min(de/e, xmax);
  G4double ds = s0*(del*G4Exp(G4Log(xmin)*bt) - bt*(xmin - 0.25*xmin*xmin));
  G4double e1 = e*(1. - xmin);
  
  //G4cout << "e1= " << e1 << G4endl;
  if(e1 < emax && s0*G4UniformRand()<ds) { 
    x = xmin*G4Exp(G4Log(G4UniformRand())/bt);
  } else {    

    x = xmin;
    G4double s1 = crossBornPerElectron->Value(e1);
    G4double w1 = bt*(del*G4Exp(G4Log(x)*btm1) - 1.0 + 0.5*x);
    G4double grej = s1*w1;
    G4double f;
    /*
     G4cout << "e(GeV)= " << e/GeV << " epeak(GeV)= " << epeak/GeV 
           << " s1= " << s1 << " w1= " << w1 
           << " grej= " << grej << G4endl;
    */
    // Above emax cross section is const
    if(e1 > emax) {
      x  = 0.5*(1. - (emax*emax)/(e*e));
      G4double s2 = crossBornPerElectron->Value(emax);
      G4double w2 = bt*(del*G4Exp(G4Log(x)*btm1) - 1.0 + 0.5*x);
      grej = s2*w2;
      //G4cout << "emax= " << emax << " s2= " << s2 << " w2= " << w2 
      // << " grej= " << grej << G4endl;
    }

    if(e1 > epeak) {
      x = 0.5*(1.0 - (epeak*epeak)/(e*e));
      G4double s2 = crossBornPerElectron->Value(epeak);
      G4double w2 = bt*(del*G4Exp(G4Log(x)*btm1) - 1.0 + 0.5*x);
      grej = std::max(grej,s2*w2);
      //G4cout << "epeak= " << epeak << " s2= " << s2 << " w2= " << w2 
      //     << " grej= " << grej << G4endl;
    }
    G4int ii = 0;
    const G4int iimax = 1000;
    do {
      x = xmin + G4UniformRand()*(xmax - xmin);
      
      G4double s2 = crossBornPerElectron->Value(sqrt(1.0 - 2*x)*e);
      G4double w2 = bt*(del*G4Exp(G4Log(x)*btm1) - 1.0 + 0.5*x);
      /*
      G4cout << "x= " << x << " xmin= " << xmin << " xmax= " << xmax
           << " s2= " << s2 << " w2= " << w2 << G4endl;
      */
      f = s2*w2;
      if(f > grej) {
	G4cout << "G4DynamicParticle* G4eeToHadronsModel:WARNING "
	       << f << " > " << grej << " majorant is`small!" 
	       << G4endl; 
      }
      if(++ii >= iimax) { break; }
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (f < grej*G4UniformRand());
  }

  G4ThreeVector dir(0.0,0.0,1.0);
  if(G4UniformRand() > 0.5) { dir.set(0.0,0.0,-1.0); }
  //G4cout << "Egamma(MeV)= " << x*e <<  " " << dir << G4endl; 
  gamma = new G4DynamicParticle(theGamma,dir,x*e);
  return gamma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
