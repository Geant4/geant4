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
// File name:     G4UrbanFluctuation
//
// Author:        V. Ivanchenko for Laszlo Urban
// 
// Creation date: 14.02.2022
//
// Modifications: 
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UrbanFluctuation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Material.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UrbanFluctuation::G4UrbanFluctuation(const G4String& nam)
 : G4UniversalFluctuation(nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UrbanFluctuation::~G4UrbanFluctuation() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UrbanFluctuation::SampleGlandz(CLHEP::HepRandomEngine* rndmEngineF,
                                          const G4Material* material,
					  const G4double tcut)
{
  if (material != lastMaterial) {
    auto ioni = material->GetIonisation();
    f1Fluct = ioni->GetF1fluct();
    f2Fluct = ioni->GetF2fluct();
    e1Fluct = ioni->GetEnergy1fluct();
    e2Fluct = ioni->GetEnergy2fluct();
    e1LogFluct = ioni->GetLogEnergy1fluct();
    e2LogFluct = ioni->GetLogEnergy2fluct();
    esmall = 0.5*std::sqrt(e0*ipotFluct);  
    lastMaterial = material;   
  }

  G4double a1(0.0), a2(0.0), a3(0.0);
  G4double loss = 0.0;
  G4double e1 = e1Fluct;
  G4double e2 = e2Fluct;

  if(tcut > ipotFluct) {
    if(w2 > ipotLogFluct)  {
      if(w2 > e2LogFluct) {
	const G4double C = meanLoss*(1.-rate)/(w2-ipotLogFluct);
	a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
	a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
      } else {
	a1 = meanLoss*(1.-rate)/e1;
      }
      if(a1 < a0) {
        const G4double fwnow = 0.5+(fw-0.5)*std::sqrt(a1/a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 *= fw;
      }
    }   
  }

  const G4double w1 = tcut/e0;
  a3 = rate*meanLoss*(tcut-e0)/(e0*tcut*G4Log(w1));
  if(a1+a2 <= 0.) { a3 /= rate; }
 
  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont  
  G4double emean = 0.;
  G4double sig2e = 0.;

  // excitation of type 1
  if(a1 > 0.0) { AddExcitation(rndmEngineF, a1, e1, emean, loss, sig2e); }

  // excitation of type 2
  if(a2 > 0.0) { AddExcitation(rndmEngineF, a2, e2, emean, loss, sig2e); }

  if(sig2e > 0.0) { SampleGauss(rndmEngineF, emean, sig2e, loss); }

  // ionisation 
  if(a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if(a3 > nmaxCont) {
      alfa = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
      const G4double alfa1  = alfa*G4Log(alfa)/(alfa-1.);
      const G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
      emean += namean*e0*alfa1;
      sig2e += e0*e0*namean*(alfa-alfa1*alfa1);
      p3 -= namean;
    }

    const G4double w3 = alfa*e0;
    if(tcut > w3) {
      const G4double w = (tcut-w3)/tcut;
      const G4int nnb = (G4int)G4Poisson(p3);
      if(nnb > 0) {
        if(nnb > sizearray) {
          sizearray = nnb;
          delete [] rndmarray;
          rndmarray = new G4double[nnb];
        }
        rndmEngineF->flatArray(nnb, rndmarray);
        for (G4int k=0; k<nnb; ++k) { loss += w3/(1.-w*rndmarray[k]); }
      }
    }
    if(sig2e > 0.0) { SampleGauss(rndmEngineF, emean, sig2e, loss); }
  }
  //G4cout << "### loss=" << loss << G4endl;
  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
