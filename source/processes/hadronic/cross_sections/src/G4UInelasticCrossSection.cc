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
// 12.08.06 V.Ivanchenko - first implementation
//
//


#include "G4UInelasticCrossSection.hh"

#include "G4ParticleTable.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4UPiNuclearCrossSection.hh"
#include "G4HadronCrossSections.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UInelasticCrossSection::G4UInelasticCrossSection(const G4ParticleDefinition* p) 
{
  fGlauber = new G4GlauberGribovCrossSection();
  fGhadProton = 0;
  if(p == G4Proton::Proton())
    fGhadProton = new G4ProtonInelasticCrossSection();
  fGhadNeutron= 0;
  if(p == G4Neutron::Neutron())
    fGhadNeutron = new G4NeutronInelasticCrossSection();
  fUPi = 0;
  if(p == G4PionPlus::PionPlus() || p == G4PionMinus::PionMinus())
    fUPi = new G4UPiNuclearCrossSection();
  fGheisha = new G4HadronCrossSections();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UInelasticCrossSection::~G4UInelasticCrossSection()
{
  delete fGlauber;
  delete fGhadProton;
  delete fGhadNeutron;
  delete fUPi;
  delete fGheisha;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4UInelasticCrossSection::IsApplicable(const G4DynamicParticle* dp, 
					      const G4Element*  elm)
{
  G4bool res = false;
  idx = 0;
  if(fGlauber->IsApplicable(dp,elm)) {
    res = true;
    idx = 1;
  } else if(fGhadProton && fGhadProton->IsApplicable(dp,elm)) {
    res = true;
    idx = 2;
  } else if(fGhadNeutron && fGhadNeutron->IsApplicable(dp,elm)) {
    res = true;
    idx = 3;
  } else if(fUPi && fUPi->IsApplicable(dp,elm)) {
    res = true;
    idx = 4;
  } else if(fGheisha->IsApplicable(dp,elm)) {
    res = true;
    idx = 5;
  }
  if(verboseLevel > 1) 
    G4cout << "G4UInelasticCrossSection::IsApplicable idx= " << idx << "  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << G4endl;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UInelasticCrossSection::GetCrossSection(const G4DynamicParticle* dp, 
						   const G4Element* elm, 
						   G4double temp)
{
  G4double cross = 0.0;
  if(idx == 1) {
    cross = fGlauber->GetCrossSection(dp,elm,temp);
    cross = fGlauber->GetInelasticGlauberGribovXsc();
  }
  else if(idx == 2) cross = fGhadProton->GetCrossSection(dp,elm,temp);
  else if(idx == 3) cross = fGhadNeutron->GetCrossSection(dp,elm,temp);
  else if(idx == 4) cross = fUPi->GetInelasticCrossSection(dp,elm);
  else if(idx == 5) cross = fGheisha->GetInelasticCrossSection(dp,elm);

  if(verboseLevel > 1) 
    G4cout << "G4UInelasticCrossSection::GetCrossSection: idx= " << idx << "  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in " << elm->GetName()
	   << G4endl;

  return cross;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4UInelasticCrossSection: uses Glauber-Gribov formula"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


