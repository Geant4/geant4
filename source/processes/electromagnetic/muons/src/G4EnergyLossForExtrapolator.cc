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
// $Id: G4EnergyLossForExtrapolator.cc 97392 2016-06-02 10:10:32Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:    G4EnergyLossForExtrapolator
//  
// Description:  This class provide calculation of energy loss, fluctuation, 
//               and msc angle
//
// Author:       09.12.04 V.Ivanchenko 
//
// Modification: 
// 08-04-05 Rename Propogator -> Extrapolator (V.Ivanchenko)
// 16-03-06 Add muon tables and fix bug in units (V.Ivanchenko)
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 12-05-06 SEt linLossLimit=0.001 (VI)
//
//----------------------------------------------------------------------------
//
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EnergyLossForExtrapolator.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TablesForExtrapolator* G4EnergyLossForExtrapolator::tables = nullptr;

G4EnergyLossForExtrapolator::G4EnergyLossForExtrapolator(G4int verb)
  : maxEnergyTransfer(DBL_MAX), verbose(verb)
{
  currentParticle = nullptr;
  currentMaterial = nullptr;

  linLossLimit = 0.01;
  emin         = 1.*MeV;
  emax         = 10.*TeV;
  nbins        = 70;

  nmat = index = 0;

  mass = charge2 = electronDensity = radLength = bg2 = beta2 
    = kineticEnergy = tmax = 0.0;
  gam = 1.0;

  idxDedxElectron = idxDedxPositron = idxDedxMuon = idxDedxProton 
    = idxRangeElectron = idxRangePositron = idxRangeMuon = idxRangeProton
    = idxInvRangeElectron = idxInvRangePositron = idxInvRangeMuon
    = idxInvRangeProton = idxMscElectron = 0;
 
  electron = positron = proton = muonPlus = muonMinus = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EnergyLossForExtrapolator::EnergyAfterStep(G4double kinEnergy, 
					     G4double stepLength, 
					     const G4Material* mat, 
					     const G4ParticleDefinition* part)
{
  if(0 == nmat) { Initialisation(); }
  G4double kinEnergyFinal = kinEnergy;
  if(SetupKinematics(part, mat, kinEnergy)) {
    G4double step = TrueStepLength(kinEnergy,stepLength,mat,part);
    G4double r  = ComputeRange(kinEnergy,part);
    if(r <= step) {
      kinEnergyFinal = 0.0;
    } else if(step < linLossLimit*r) {
      kinEnergyFinal -= step*ComputeDEDX(kinEnergy,part);
    } else {  
      G4double r1 = r - step;
      kinEnergyFinal = ComputeEnergy(r1,part);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EnergyLossForExtrapolator::EnergyBeforeStep(G4double kinEnergy, 
					      G4double stepLength, 
					      const G4Material* mat, 
					      const G4ParticleDefinition* part)
{
  //G4cout << "G4EnergyLossForExtrapolator::EnergyBeforeStep" << G4endl;
  if(0 == nmat) { Initialisation(); }
  G4double kinEnergyFinal = kinEnergy;

  if(SetupKinematics(part, mat, kinEnergy)) {
    G4double step = TrueStepLength(kinEnergy,stepLength,mat,part);
    G4double r  = ComputeRange(kinEnergy,part);

    if(step < linLossLimit*r) {
      kinEnergyFinal += step*ComputeDEDX(kinEnergy,part);
    } else {  
      G4double r1 = r + step;
      kinEnergyFinal = ComputeEnergy(r1,part);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EnergyLossForExtrapolator::TrueStepLength(G4double kinEnergy, 
					    G4double stepLength,
					    const G4Material* mat, 
					    const G4ParticleDefinition* part)
{
  G4double res = stepLength;
  if(0 == nmat) { Initialisation(); }
  //G4cout << "## G4EnergyLossForExtrapolator::TrueStepLength L= " << res 
  //	 <<  "  " << part->GetParticleName() << G4endl;
  if(SetupKinematics(part, mat, kinEnergy)) {
    if(part == electron || part == positron) {
      G4double x = stepLength*
	ComputeValue(kinEnergy, GetPhysicsTable(fMscElectron), idxMscElectron);
      //G4cout << " x= " << x << G4endl;
      if(x < 0.2)         { res *= (1.0 + 0.5*x + x*x/3.0); }
      else if(x < 0.9999) { res = -G4Log(1.0 - x)*stepLength/x; }
      else                { res = ComputeRange(kinEnergy,part); }
    } else {
      res = ComputeTrueStep(mat,part,kinEnergy,stepLength);
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool 
G4EnergyLossForExtrapolator::SetupKinematics(const G4ParticleDefinition* part, 
					     const G4Material* mat, 
					     G4double kinEnergy)
{
  if(0 == nmat) { Initialisation(); }
  if(!part || !mat || kinEnergy < keV) { return false; }
  G4bool flag = false;
  if(part != currentParticle) {
    flag = true;
    currentParticle = part;
    mass = part->GetPDGMass();
    G4double q = part->GetPDGCharge()/eplus;
    charge2 = q*q;
  }
  if(mat != currentMaterial) {
    G4int i = mat->GetIndex();
    if(i >= nmat) {
      G4cout << "### G4EnergyLossForExtrapolator WARNING:index i= " 
	     << i << " is out of table - NO extrapolation" << G4endl;
    } else {
      flag = true;
      currentMaterial = mat;
      electronDensity = mat->GetElectronDensity();
      radLength       = mat->GetRadlen();
      index           = i;
    }
  }
  if(flag || kinEnergy != kineticEnergy) {
    kineticEnergy = kinEnergy;
    G4double tau  = kinEnergy/mass;

    gam   = tau + 1.0;
    bg2   = tau * (tau + 2.0);
    beta2 = bg2/(gam*gam);
    tmax  = kinEnergy;
    if(part == electron) tmax *= 0.5;
    else if(part != positron) {
      G4double r = electron_mass_c2/mass;
      tmax = 2.0*bg2*electron_mass_c2/(1.0 + 2.0*gam*r + r*r);
    }
    if(tmax > maxEnergyTransfer) { tmax = maxEnergyTransfer; }
  }
  return true;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* 
G4EnergyLossForExtrapolator::FindParticle(const G4String& name)
{
  const G4ParticleDefinition* p = nullptr;
  if(name != currentParticleName) {
    p = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if(!p) {
      G4cout << "### G4EnergyLossForExtrapolator WARNING: "
	     << "FindParticle() fails to find " 
	     << name << G4endl;
    }
  } else {
    p = currentParticle;
  }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EnergyLossForExtrapolator::ComputeDEDX(G4double ekin, 
					 const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)  { 
    x = ComputeValue(ekin, GetPhysicsTable(fDedxElectron), idxDedxElectron);
  } else if(part == positron) {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxPositron), idxDedxPositron);
  } else if(part == muonPlus || part == muonMinus) {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxMuon), idxDedxMuon);
  } else {
    G4double e = ekin*proton_mass_c2/mass;
    x = ComputeValue(e, GetPhysicsTable(fDedxProton), idxDedxProton)*charge2;
  }
  return x;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EnergyLossForExtrapolator::ComputeRange(G4double ekin, 
					  const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron) { 
    x = ComputeValue(ekin, GetPhysicsTable(fRangeElectron), idxRangeElectron);
  } else if(part == positron) {
    x = ComputeValue(ekin, GetPhysicsTable(fRangePositron), idxRangePositron);
  } else if(part == muonPlus || part == muonMinus) { 
    x = ComputeValue(ekin, GetPhysicsTable(fRangeMuon), idxRangeMuon);
  } else {
    G4double massratio = proton_mass_c2/mass;
    G4double e = ekin*massratio;
    x = ComputeValue(e, GetPhysicsTable(fRangeProton), idxRangeProton)
      /(charge2*massratio);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EnergyLossForExtrapolator::ComputeEnergy(G4double range, 
					   const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron) {
    x = ComputeValue(range, GetPhysicsTable(fInvRangeElectron), 
		     idxInvRangeElectron);
  } else if(part == positron) {
    x = ComputeValue(range, GetPhysicsTable(fInvRangePositron), 
		     idxInvRangePositron);
  } else if(part == muonPlus || part == muonMinus) {
    x = ComputeValue(range, GetPhysicsTable(fInvRangeMuon), idxInvRangeMuon);
  } else {
    G4double massratio = proton_mass_c2/mass;
    G4double r = range*massratio*charge2;
    x = ComputeValue(r, GetPhysicsTable(fInvRangeProton), 
		     idxInvRangeProton)/massratio;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::Initialisation()
{
  if(verbose>1) {
    G4cout << "### G4EnergyLossForExtrapolator::Initialisation" << G4endl;
  }
  currentParticle = nullptr;
  currentMaterial = nullptr;
  kineticEnergy   = 0.0;

  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton   = G4Proton::Proton();
  muonPlus = G4MuonPlus::MuonPlus();
  muonMinus= G4MuonMinus::MuonMinus();

  nmat = G4Material::GetNumberOfMaterials();
  currentParticleName = "";
  BuildTables(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4AutoLock.hh"
namespace { G4Mutex G4EnergyLossForExtrapolatorMutex = G4MUTEX_INITIALIZER; }

void G4EnergyLossForExtrapolator::BuildTables()
{
  G4AutoLock l(&G4EnergyLossForExtrapolatorMutex);

  if(!tables) {
    if(verbose > 0) {
      G4cout << "### G4EnergyLossForExtrapolator::BuildTables for "
	     << nmat << " materials Nbins= " << nbins 
	     << " Emin(MeV)= " << emin << "  Emax(MeV)= " << emax << G4endl;
    }
    tables = new G4TablesForExtrapolator(verbose, nbins, emin, emax);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolator::~G4EnergyLossForExtrapolator()
{
  G4AutoLock l(&G4EnergyLossForExtrapolatorMutex);
  delete tables;
  tables = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

