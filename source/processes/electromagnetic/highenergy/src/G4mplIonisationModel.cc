//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4mplIonisationModel.cc,v 1.1 2006-10-25 17:37:44 vnivanch Exp $
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
//
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4mplIonisationModel.hh"
#include "Randomize.hh"
#include "G4Monopole.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForLoss.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4mplIonisationModel::G4mplIonisationModel(const G4ParticleDefinition*, const G4String& nam)
  : G4VEmModel(nam),
  twoln10(2.0*log(10.0)),
  beta2low(0.0001),
  beta2lim(0.01),
  bg2lim(beta2lim*(1.0 + beta2lim))
{
  monopole     = G4Monopole::Monopole();
  mass         = monopole->GetPDGMass();
  magCharge    = monopole->MagneticCharge()/eplus;
  nmpl         = G4int(abs(magCharge)/68.0);
  if(nmpl > 6)      nmpl = 6;
  else if(nmpl < 1) nmpl = 1;
  G4double x   = 45.0*GeV*G4double(nmpl)/cm;
  factlow      = x*x;
  chargeSquare = magCharge*magCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4mplIonisationModel::~G4mplIonisationModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4mplIonisationModel::Initialise(const G4ParticleDefinition*,
				      const G4DataVector&)
{
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
  G4double tau   = kineticEnergy/mass;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);

  G4double dedx0 = factlow*abs(beta2);

  if(beta2 > beta2low) {

    G4double b2 = beta2;
    if(beta2 < beta2lim) {
      beta2= beta2lim;
      bg2  = bg2lim;
    }

    G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
    G4double cden  = material->GetIonisation()->GetCdensity();
    G4double mden  = material->GetIonisation()->GetMdensity();
    G4double aden  = material->GetIonisation()->GetAdensity();
    G4double x0den = material->GetIonisation()->GetX0density();
    G4double x1den = material->GetIonisation()->GetX1density();

    G4double eDensity = material->GetElectronDensity();

    G4double dedx = 2.0*log(2.0*electron_mass_c2*bg2/eexc) - 1.0;

    G4double  k = 0.406;
    if(nmpl > 1) k = 0.346;
    const G4double B[7] = { 0.0, 0.248, 0.672, 1.022,  1.243, 1.464,  1.685}; 

    dedx += k - B[nmpl];

    // density correction
    G4double x = log(bg2)/twoln10;
    if ( x >= x0den ) {
      dedx -= twoln10*x - cden ;
      if ( x < x1den ) dedx -= aden*pow((x1den-x),mden) ;
    }

    // now compute the total ionization loss

    if (dedx < 0.0) dedx = 0.0 ;

    dedx *= twopi_mc2_rcl2*chargeSquare*eDensity;

    // extrapolate between two formula
    if(beta2 < beta2lim) {
      x = log(dedx0) + log(dedx/dedx0)*log(b2/beta2low)/log(beta2lim/beta2low);
      dedx = exp(x);
    }
    dedx0 = dedx;
  }

  return dedx0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
