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
// Geant4 class G4HadXSHelper
//
// Author V.Ivanchenko 18.05.2022
//

#include "G4HadXSHelper.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static const G4double scale = 10./G4Log(10.);

std::vector<G4double>* 
G4HadXSHelper::FindCrossSectionMax(G4HadronicProcess* p,
                                   const G4ParticleDefinition* part,
                                   const G4double emin,
                                   const G4double emax)
{
  std::vector<G4double>* ptr = nullptr;
  if(nullptr == p || nullptr == part) { return ptr; }
  /*
  G4cout << "G4HadXSHelper::FindCrossSectionMax for "
	 << p->GetProcessName() << " and " << part->GetParticleName() << G4endl;
  */

  const G4MaterialTable* theMatTable = G4Material::GetMaterialTable();
  size_t n = G4Material::GetNumberOfMaterials();
  ptr = new std::vector<G4double>;
  ptr->resize(n, DBL_MAX);

  G4bool isPeak = false;
  G4double ee = G4Log(emax/emin);
  G4int nbin = G4lrint(ee*scale);
  if(nbin < 4) { nbin = 4; }
  G4double x = G4Exp(ee/nbin);

  // first loop on existing vectors
  for (size_t i=0; i<n; ++i) {
    auto mat = (*theMatTable)[i];
    G4double sm = 0.0;
    G4double em = 0.0;
    G4double e = emin;
    for(G4int j=0; j<=nbin; ++j) {
      G4double sig = p->ComputeCrossSection(part, mat, e);
      if(sig >= sm) {
	em = e;
	sm = sig;
	e = (j+1 < nbin) ? e*x : emax;
      } else {
	isPeak = true;
	(*ptr)[i] = em;
	break;
      }
    }
    //G4cout << i << ".  em=" << em << " sm=" << sm << G4endl;
  }
  // there is no peak for any couple
  if(!isPeak) {
    delete ptr;
    ptr = nullptr;
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<G4TwoPeaksHadXS*>* 
G4HadXSHelper::FillPeaksStructure(G4HadronicProcess* p,
				  const G4ParticleDefinition* part,
				  const G4double emin,
				  const G4double emax)
{
  std::vector<G4TwoPeaksHadXS*>* ptr = nullptr;
  if(nullptr == p) { return ptr; }

  /*
  G4cout << "G4HadXSHelper::FillPeaksStructure for "
	 << p->GetProcessName() << " and " << part->GetParticleName() << G4endl;
  */
  const G4MaterialTable* theMatTable = G4Material::GetMaterialTable();
  size_t n = G4Material::GetNumberOfMaterials();
  ptr = new std::vector<G4TwoPeaksHadXS*>;
  ptr->resize(n, nullptr);

  G4double e1peak, e1deep, e2peak, e2deep, e3peak;
  G4bool isDeep = false;

  G4double ee = G4Log(emax/emin);
  G4int nbin = G4lrint(ee*scale);
  if(nbin < 4) { nbin = 4; }
  G4double fact = G4Exp(ee/nbin);

  for (size_t i=0; i<n; ++i) {
    auto mat = (*theMatTable)[i];
    G4double e = emin/fact;

    G4double xs = 0.0;
    e1peak = e1deep = e2peak = e2deep = e3peak = DBL_MAX;
    for(G4int j=0; j<=nbin; ++j) {
      e = (j+1 < nbin) ? e*fact : emax;
      G4double ss = p->ComputeCrossSection(part, mat, e);
      // find out 1st peak
      if(e1peak == DBL_MAX) {
	if(ss >= xs) {
	  xs = ss;
	  ee = e;
	  continue;
	} else {
	  e1peak = ee;
	}
      }
      // find out the deep
      if(e1deep == DBL_MAX) {
	if(ss <= xs) {
	  xs = ss;
	  ee = e;
	  continue;
	} else {
	  e1deep = ee;
	  isDeep = true;
	}
      }
      // find out 2nd peak
      if(e2peak == DBL_MAX) {
	if(ss >= xs) {
	  xs = ss;
	  ee = e;
	  continue;
	} else {
	  e2peak = ee;
	}
      }
      if(e2deep == DBL_MAX) {
	if(ss <= xs) {
	  xs = ss;
	  ee = e;
	  continue;
	} else {
	  e2deep = ee;
	  break;
	}
      }
      // find out 3d peak
      if(e3peak == DBL_MAX) {
	if(ss >= xs) {
	  xs = ss;
	  ee = e;
	  continue;
	} else {
	  e3peak = ee;
	}
      }
    }
    G4TwoPeaksHadXS* x = (*ptr)[i];
    if(nullptr == x) { 
      x = new G4TwoPeaksHadXS(); 
      (*ptr)[i] = x;
    }
    x->e1peak = e1peak;
    x->e1deep = e1deep;
    x->e2peak = e2peak;
    x->e2deep = e2deep;
    x->e3peak = e3peak;
    //G4cout << "coupleIdx=" << i << " " << e1peak << " " << e1deep
    //  << " " << e2peak << " " << e2deep << " " << e3peak << G4endl;
  }
  // case of no 1st peak in all vectors
  if(!isDeep) {
    delete ptr;
    ptr = nullptr;
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
