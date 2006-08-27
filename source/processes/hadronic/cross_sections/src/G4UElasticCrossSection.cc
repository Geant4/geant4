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


#include "G4UElasticCrossSection.hh"

#include "G4ParticleTable.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4UPiNuclearCrossSection.hh"
#include "G4HadronCrossSections.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UElasticCrossSection::G4UElasticCrossSection() 
{
  fGlauber = new G4GlauberGribovCrossSection();
  fUPi     = new G4UPiNuclearCrossSection();
  fGheisha = new G4HadronCrossSections();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UElasticCrossSection::~G4UElasticCrossSection()
{
  delete fGlauber;
  delete fGheisha;
  delete fUPi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4UElasticCrossSection::IsApplicable(const G4DynamicParticle* dp, 
					    const G4Element*  elm)
{
  G4bool res = false;
  idx = 0;
  if(fGlauber->IsApplicable(dp,elm)) {
    res = true;
    idx = 1;
  } else if(fUPi->IsApplicable(dp,elm)) {
    res = true;
    idx = 2;
  } else if(fGheisha->IsApplicable(dp,elm)) {
    res = true;
    idx = 3;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UElasticCrossSection::GetCrossSection(const G4DynamicParticle* dp, 
						 const G4Element* elm, 
						 G4double temp)
{
  G4double cross = 0.0;
  if(idx == 1) {
    cross = fGlauber->GetCrossSection(dp,elm,temp);
    cross = fGlauber->GetElasticGlauberGribovXsc();
  } 
  else if(idx == 2) cross = fUPi->GetElasticCrossSection(dp,elm);
  else if(idx == 3) cross = fGheisha->GetElasticCrossSection(dp,elm);

  return cross;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UElasticCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UElasticCrossSection::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4UElasticCrossSection: uses Glauber-Gribov formula"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


