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
// $Id: G4BGGPionElasticXS.cc,v 1.1 2007/03/13 15:19:30 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BGGPionElasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 01.10.2003
// Modifications:
//
//
// -------------------------------------------------------------------
//

#include "G4BGGPionElasticXS.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4UPiNuclearCrossSection.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGPionElasticXS::G4BGGPionElasticXS(const G4ParticleDefinition* p) 
{
  verboseLevel = 0;
  thEnergy     = 100.*GeV;
  if(p == G4PionPlus::PionPlus() || p == G4PionMinus::PionMinus()) {
    fPion = new G4UPiNuclearCrossSection();
    fGlauber = new G4GlauberGribovCrossSection();
    particle = p;
    Initialise();
  } else {
    fPion = 0;
    fGlauber = 0;
    particle = 0;
    if(p) G4cout << "### G4BGGPionElasticXS WARNING: is not applicable to " 
		 << p->GetParticleName()
		 << G4endl;
    else  G4cout << "### G4BGGPionElasticXS WARNING: particle is not defined " 
		 << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGPionElasticXS::~G4BGGPionElasticXS()
{
  delete fGlauber;
  delete fPion;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGPionElasticXS::GetIsoZACrossSection(const G4DynamicParticle* dp, 
						    G4double Z,
						    G4double A, 
						    G4double)
{
  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int iz = G4int(Z + 0.5);
  if(iz > 92) iz = 92;

  if(ekin > thEnergy) {
    cross = theFac[iz]*fGlauber->GetElasticGlauberGribov(dp, Z, A);
  } else {
    cross = fPion->GetElasticCrossSection(dp, Z, A);
  }

  if(verboseLevel > 1) 
    G4cout << "G4BGGPionElasticXS::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGPionElasticXS::BuildPhysicsTable(const G4ParticleDefinition&)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGPionElasticXS::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4BGGPionElasticXS:"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGPionElasticXS::Initialise() 
{
  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(particle);
  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(part, mom, thEnergy);

  G4NistManager* nist = G4NistManager::Instance();
  G4double A = nist->GetAtomicMassAmu(2);

  G4double csup, csdn;

  if(verboseLevel > 0) G4cout << "### G4BGGPionElasticXS::Initialise for "
			      << particle->GetParticleName() << G4endl;

  for(G4int iz=2; iz<93; iz++) {

    G4double Z = G4double(iz);
    A = nist->GetAtomicMassAmu(iz);

    csup = fGlauber->GetElasticGlauberGribov(&dp, Z, A);
    csdn = fPion->GetElasticCrossSection(&dp, Z, A);

    theFac[iz] = csdn/csup;
    if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				<< " factor= " << theFac[iz] << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


