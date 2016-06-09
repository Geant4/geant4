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
// $Id: G4PEEffectModel.cc,v 1.3 2005/08/11 10:03:33 maire Exp $
// GEANT4 tag $Name: geant4-07-01-patch-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PEEffectModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 21.03.2005
//
// Modifications:
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PEEffectModel.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4PEEffectModel::G4PEEffectModel(const G4ParticleDefinition*,
					 const G4String& nam)
  : G4VEmModel(nam),isInitialized(false)
{
  theGamma    = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PEEffectModel::~G4PEEffectModel()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PEEffectModel::Initialise(const G4ParticleDefinition*,
				 const G4DataVector&)
{
  if(isInitialized) return;
  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForGamma();

  fminimalEnergy = 1.0*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4PEEffectModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* aDynamicPhoton,
                                   G4double,
                                   G4double)
{
  const G4Material* aMaterial = couple->GetMaterial();

  G4double energy = aDynamicPhoton->GetKineticEnergy();
  G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();

  // select randomly one element constituing the material.
  const G4Element* anElement = SelectRandomAtom(aMaterial,theGamma,energy);
  
  //
  // Photo electron
  //
  std::vector<G4DynamicParticle*>* fvect = new std::vector<G4DynamicParticle*>;

  // Select atomic shell
  G4int nShells = anElement->GetNbOfAtomicShells();
  G4int i  = 0;  
  while ((i<nShells) && (energy<anElement->GetAtomicShell(i))) i++;

  // no shell available
  if (i == nShells) return fvect;
  
  G4double bindingEnergy  = anElement->GetAtomicShell(i);
  G4double ElecKineEnergy = energy - bindingEnergy;

  if (ElecKineEnergy > fminimalEnergy)
    {
     // direction of the photo electron
     //
     G4double cosTeta = ElecCosThetaDistribution(ElecKineEnergy);
     G4double sinTeta = sqrt(1.-cosTeta*cosTeta);
     G4double Phi     = twopi * G4UniformRand();
     G4double dirx = sinTeta*cos(Phi),diry = sinTeta*sin(Phi),dirz = cosTeta;
     G4ThreeVector ElecDirection(dirx,diry,dirz);
     ElecDirection.rotateUz(PhotonDirection);
     //
     G4DynamicParticle* aParticle = new G4DynamicParticle (
                       theElectron,ElecDirection, ElecKineEnergy);
     fvect->push_back(aParticle);
    }

  fParticleChange->ProposeTrackStatus(fStopAndKill);   
  fParticleChange->ProposeLocalEnergyDeposit(bindingEnergy);
  return fvect;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PEEffectModel::ElecCosThetaDistribution(G4double kineEnergy)
{
  // Compute Theta distribution of the emitted electron, with respect to the
  // incident Gamma.
  // The Sauter-Gavrila distribution for the K-shell is used.
  //
  G4double costeta = 1.;
  G4double gamma   = 1. + kineEnergy/electron_mass_c2;
  if (gamma > 5.) return costeta;
  G4double beta  = sqrt(gamma*gamma-1.)/gamma;
  G4double b     = 0.5*gamma*(gamma-1.)*(gamma-2);

  G4double rndm,term,greject,grejsup;
  if (gamma < 2.) grejsup = gamma*gamma*(1.+b-beta*b);
  else            grejsup = gamma*gamma*(1.+b+beta*b);

  do { rndm = 1.-2*G4UniformRand();
       costeta = (rndm+beta)/(rndm*beta+1.);
       term = 1.-beta*costeta;
       greject = (1.-costeta*costeta)*(1.+b*term)/(term*term);
  } while(greject < G4UniformRand()*grejsup);

  return costeta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
