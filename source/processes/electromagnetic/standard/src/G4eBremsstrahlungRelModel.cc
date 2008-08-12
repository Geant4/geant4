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
// $Id: G4eBremsstrahlungRelModel.cc,v 1.1 2008-08-12 17:51:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eBremsstrahlungRelModel
//
// Author:        Andreas Schaelicke 
//
// Creation date: 12.08.2008
//
// Modifications:
//
// Main References:
//  S.Klein,  Rev. Mod. Phys. 71 (1999) 1501.
//  T.Stanev et.al., Phys. Rev. D25 (1982) 1291.
//  M.L.Ter-Mikaelian, High-energy Electromagnetic Processes in Condensed Media, Wiley, 1972.
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eBremsstrahlungRelModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::xgi[]={ 0.0199, 0.1017, 0.2372, 0.4083,
					    0.5917, 0.7628, 0.8983, 0.9801 };
G4double G4eBremsstrahlungRelModel::wgi[]={ 0.0506, 0.1112, 0.1569, 0.1813,
					    0.1813, 0.1569, 0.1112, 0.0506 };

using namespace std;

G4eBremsstrahlungRelModel::G4eBremsstrahlungRelModel(const G4ParticleDefinition* p,
						     const G4String& nam)
  : G4VEmModel(nam),
    particle(0),
    isElectron(true),
    //    MigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length/pi),
    MigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length*4.0*pi),
    LPMconstant(fine_structure_const*electron_mass_c2*electron_mass_c2/(4.*pi*hbarc)),
    bremFactor(fine_structure_const*classic_electr_radius*classic_electr_radius*4./3.),
    isInitialised(false)
{
  if(p) SetParticle(p);
  theGamma = G4Gamma::Gamma();
  minThreshold = 1.0*keV;
  highEnergyTh = DBL_MAX;
  hydrogenEnergyTh = 0.0;
  nist = G4NistManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlungRelModel::~G4eBremsstrahlungRelModel()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungRelModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  particleMass = p->GetPDGMass();
  if(p == G4Electron::Electron()) isElectron = true;
  else                            isElectron = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::MinEnergyCut(const G4ParticleDefinition*,
						 const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungRelModel::SetupForMaterial(const G4ParticleDefinition*,
						 const G4Material* mat)
{
  densityFactor = mat->GetElectronDensity()*MigdalConstant;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungRelModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector& cuts)
{
  // *** update flags from losstablemanager *** 
  G4LossTableManager* man = G4LossTableManager::Instance(); 
  SetEnergyThreshold(man->BremsstrahlungTh());

  if(p) SetParticle(p);

  highKinEnergy = HighEnergyLimit();
  lowKinEnergy  = LowEnergyLimit();

  currentZ = 0;

  InitialiseElementSelectors(p, cuts);

  if(isInitialised) return;

  if(pParticleChange) {
    fParticleChange = reinterpret_cast<G4ParticleChangeForLoss*>(pParticleChange);
  } else {
    fParticleChange = new G4ParticleChangeForLoss();
  }
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::ComputeDEDXPerVolume(
					     const G4Material* material,
                                             const G4ParticleDefinition* p,
                                                   G4double kineticEnergy,
                                                   G4double cutEnergy)
{
  if(!particle) SetParticle(p);
  if(kineticEnergy < lowKinEnergy) return 0.0;
  G4double cut = std::min(cutEnergy, kineticEnergy);
  if(cut == 0.0) return 0.0;

  SetupForMaterial(particle, material);

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();

  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;

  G4double dedx = 0.0;

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    SetCurrentElement((*theElementVector)[i]->GetZ());
    dedx += theAtomicNumDensityVector[i]*currentZ*currentZ*ComputeBremLoss(cut);
  }
  dedx *= bremFactor;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::ComputeBremLoss(G4double cut)
{
  G4double loss = 0.0;

  // number of intervals and integration step 
  G4double vcut = cut/totalEnergy;
  G4int n = (G4int)(20*vcut) + 3;
  G4double delta = vcut/G4double(n);

  G4double e0 = 0.0;
  G4double xs; 

  // integration
  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = (e0 + xgi[i]*delta)*totalEnergy;

      if(totalEnergy > hydrogenEnergyTh*currentZ) {
	xs = ComputeRelDXSectionPerAtom(eg);
      } else {
	xs = ComputeDXSectionPerAtom(eg);
      }
      loss += eg*wgi[i]*xs;
    }
    e0 += delta;
  }

  loss *= delta*totalEnergy;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::ComputeCrossSectionPerAtom(
                                              const G4ParticleDefinition* p,
					      G4double kineticEnergy, 
					      G4double Z,   G4double,
					      G4double cutEnergy, 
					      G4double maxEnergy)
{
  if(!particle) SetParticle(p);
  if(kineticEnergy < lowKinEnergy) return 0.0;
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  G4double tmax = std::min(maxEnergy, kineticEnergy);
  if(cut >= tmax) return 0.0;

  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;
  SetCurrentElement(Z);

  G4double cross = ComputeXSectionPerAtom(cut);

  // allow partial integration
  if(tmax < kinEnergy) cross -= ComputeXSectionPerAtom(tmax);
  
  cross *= Z*Z*bremFactor;
  return cross;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4eBremsstrahlungRelModel::ComputeXSectionPerAtom(G4double cut)
{
  G4double cross = 0.0;

  // number of intervals and integration step 
  G4double vcut = log(cut/totalEnergy);
  G4double vmax = log(kinEnergy/totalEnergy);
  G4int n = (G4int)(0.45*(vmax - vcut)) + 4;
  G4double delta = (vmax - vcut)/G4double(n);

  G4double e0 = vcut;
  G4double xs; 

  // integration
  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = (e0 + xgi[i]*delta)*totalEnergy;

      if(totalEnergy > hydrogenEnergyTh*currentZ) {
	xs = ComputeRelDXSectionPerAtom(eg);
      } else {
	xs = ComputeDXSectionPerAtom(eg);
      }
      cross += eg*wgi[i]*xs;
    }
    e0 += delta;
  }

  cross *= delta;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::ComputeRelDXSectionPerAtom(G4double gammaEnergy)
  // Ultra relativistic model
{
  G4double cross = 0.0;
  if(gammaEnergy < 0.0) cross = 0.0;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungRelModel::ComputeDXSectionPerAtom(G4double gammaEnergy)
  // Relativistic model
{
  G4double cross = 0.0;
  if(gammaEnergy < 0.0) cross = 0.0;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungRelModel::SampleSecondaries(
				      std::vector<G4DynamicParticle*>* vdp, 
				      const G4MaterialCutsCouple* couple,
				      const G4DynamicParticle* dp,
				      G4double cutEnergy,
				      G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  if(kineticEnergy < lowKinEnergy) return;
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  G4double emax = std::min(maxEnergy, kineticEnergy);
  if(cut >= emax) return;

  SetupForMaterial(particle, couple->GetMaterial());

  const G4Element* elm = 
    SelectRandomAtom(couple,particle,kineticEnergy,cut,emax);
  SetCurrentElement(elm->GetZ());

  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;
  G4ThreeVector direction = dp->GetMomentumDirection();

  G4double fmax= 1.0;
  G4bool highe = true;
  if(totalEnergy < hydrogenEnergyTh*currentZ) highe = false;
 
  G4double xmin = log(cut);
  G4double xmax = log(emax);
  G4double gammaEnergy, f; 

  do {
    gammaEnergy = exp(xmin + G4UniformRand()*(xmax - xmin));
    if(highe) f = ComputeRelDXSectionPerAtom(gammaEnergy);
    else      f = ComputeDXSectionPerAtom(gammaEnergy);

    if ( f > fmax ) {
      G4cout << "### G4eBremsstrahlungRelModel Warning: Majoranta exceeded! "
	     << f << " > " << fmax
	     << " Egamma(MeV)= " << gammaEnergy
	     << " E(mEV)= " << kineticEnergy
	     << G4endl;
    }

  } while (f < fmax*G4UniformRand());

  //
  //  angles of the emitted gamma. ( Z - axis along the parent particle)
  //
  //  universal distribution suggested by L. Urban 
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > G4UniformRand()) u = - log(G4UniformRand()*G4UniformRand())/a1;
     else                          u = - log(G4UniformRand()*G4UniformRand())/a2;

  G4double theta = u*particleMass/totalEnergy;
  G4double sint = sin(theta);
  G4double phi = twopi * G4UniformRand();
  G4ThreeVector gammaDirection(sint*cos(phi),sint*sin(phi), cos(theta));
  gammaDirection.rotateUz(direction);

  // create G4DynamicParticle object for the Gamma
  G4DynamicParticle* g = new G4DynamicParticle(theGamma,gammaDirection,
                                                        gammaEnergy);
  vdp->push_back(g);
  
  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
  G4ThreeVector dir = totMomentum*direction - gammaEnergy*gammaDirection;
  direction = dir.unit();

  // energy of primary
  G4double finalE = kineticEnergy - gammaEnergy;

  // stop tracking and create new secondary instead of primary
  if(gammaEnergy > highEnergyTh) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    G4DynamicParticle* el = 
      new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
			    direction, finalE);
    vdp->push_back(el);

    // continue tracking
  } else {
    fParticleChange->SetProposedMomentumDirection(direction);
    fParticleChange->SetProposedKineticEnergy(finalE);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


