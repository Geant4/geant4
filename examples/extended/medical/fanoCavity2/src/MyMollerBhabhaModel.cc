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
/// \file medical/fanoCavity2/src/MyMollerBhabhaModel.cc
/// \brief Implementation of the MyMollerBhabhaModel class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MyMollerBhabhaModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

MyMollerBhabhaModel::MyMollerBhabhaModel(const G4ParticleDefinition* p,
                                         const G4String& nam)
  : G4MollerBhabhaModel(p,nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyMollerBhabhaModel::~MyMollerBhabhaModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MyMollerBhabhaModel::ComputeDEDXPerVolume(
                                          const G4Material* material,
                                          const G4ParticleDefinition* p,
                                                G4double kineticEnergy,
                                                G4double cutEnergy)
{
  if(!particle) SetParticle(p);
  // calculate the dE/dx due to the ionization by Seltzer-Berger formula
  
  G4double electronDensity = material->GetElectronDensity();
  G4double Zeff  = electronDensity/material->GetTotNbOfAtomsPerVolume();
  G4double th    = 0.25*sqrt(Zeff)*keV;
  G4double tkin  = kineticEnergy;
  G4bool   lowEnergy = false;
  if (kineticEnergy < th) {
    tkin = th;
    lowEnergy = true;
  }
  G4double tau   = tkin/electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double beta2 = 1. - 1./gamma2;
  //G4double bg2   = beta2*gamma2;

  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  eexc          /= electron_mass_c2;
  G4double eexc2 = eexc*eexc; 

  G4double d = min(cutEnergy, MaxSecondaryEnergy(p, tkin))/electron_mass_c2;
  G4double dedx;

  // electron
  if (isElectron) {

    dedx = log(2.0*(tau + 2.0)/eexc2) - 1.0 - beta2
         + log((tau-d)*d) + tau/(tau-d)
         + (0.5*d*d + (2.0*tau + 1.)*log(1. - d/tau))/gamma2;
   
  //positron
  } else {

    G4double d2 = d*d*0.5;
    G4double d3 = d2*d/1.5;
    G4double d4 = d3*d*3.75;
    G4double y  = 1.0/(1.0 + gam);
    dedx = log(2.0*(tau + 2.0)/eexc2) + log(tau*d)
         - beta2*(tau + 2.0*d - y*(3.0*d2 
         + y*(d - d3 + y*(d2 - tau*d3 + d4))))/tau;
  } 

  //do not apply density correction 
  //G4double cden  = material->GetIonisation()->GetCdensity();
  //G4double mden  = material->GetIonisation()->GetMdensity();
  //G4double aden  = material->GetIonisation()->GetAdensity();
  //G4double x0den = material->GetIonisation()->GetX0density();
  //G4double x1den = material->GetIonisation()->GetX1density();
  //G4double x     = log(bg2)/twoln10;
  
  //if (x >= x0den) {
  //  dedx -= twoln10*x - cden;
  //  if (x < x1den) dedx -= aden*pow(x1den-x, mden);
  //} 

  // now you can compute the total ionization loss
  dedx *= twopi_mc2_rcl2*electronDensity/beta2;
  if (dedx < 0.0) dedx = 0.0;

  // lowenergy extrapolation

  if (lowEnergy) {

    if (kineticEnergy >= lowLimit) dedx *= sqrt(tkin/kineticEnergy);
    else                           dedx *= sqrt(tkin*kineticEnergy)/lowLimit;

  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
