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
// $Id: G4LivermoreIonisationCrossSection.cc 66241 2012-12-13 18:34:42Z gunter $
//
// Author: Vladimir Ivanchenko
//
// History:
// --------
// 31 May 2011 V.Ivanchenko  Created  
// 09 Mar 2012 L.Pandola     update methods
// 
//

#include "G4LivermoreIonisationCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4eCrossSectionHandler.hh"
#include "G4SemiLogInterpolation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4LivermoreIonisationCrossSection::G4LivermoreIonisationCrossSection(
  const G4String& nam) : G4VhShellCrossSection(nam), crossSectionHandler(0)
{
  fLowEnergyLimit  = 10.0*eV;
  fHighEnergyLimit = 100.0*GeV;

  transitionManager = G4AtomicTransitionManager::Instance();

  verboseLevel = 0;
  
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreIonisationCrossSection::~G4LivermoreIonisationCrossSection()
{
  delete crossSectionHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationCrossSection::Initialise()
{
  const G4int binForFluo = 20;
  G4int nbin = G4int(std::log10(fHighEnergyLimit/fLowEnergyLimit) + 0.5);
  if(nbin <= 0) { nbin = 1; }
  nbin *= binForFluo;

  // Data on shell ionisation x-sections
  if (crossSectionHandler) { 
    crossSectionHandler->Clear();
    delete crossSectionHandler; 
  }

  G4VDataSetAlgorithm* inter = new G4SemiLogInterpolation();
  crossSectionHandler = 
    new G4eCrossSectionHandler(inter,fLowEnergyLimit,fHighEnergyLimit,nbin);
  crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");
  //G4cout << "!!!  G4LivermoreIonisationCrossSection::Initialise()" << G4endl;
}

G4double 
G4LivermoreIonisationCrossSection::CrossSection(G4int Z, G4AtomicShellEnumerator shell,
						G4double kinEnergy, G4double,
						const G4Material*) 
{
  G4double cross = 0.0;
  G4int n = G4int(shell);
  G4int nmax = std::min(9,transitionManager->NumberOfShells(Z));
  if(Z > 6 && Z < 93 && n < nmax && 
     kinEnergy >= fLowEnergyLimit && kinEnergy <= fHighEnergyLimit) {
    //G4cout << "Z= " << Z << "  n= " << n << " E(MeV)= " << kinEnergy/MeV << G4endl;
    cross = crossSectionHandler->FindValue(Z, kinEnergy, n);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4double> 
G4LivermoreIonisationCrossSection::GetCrossSection(G4int Z,
						   G4double kinEnergy,
						   G4double, G4double, 
						   const G4Material*)
{
  G4int nmax = std::min(9,transitionManager->NumberOfShells(Z));
  std::vector<G4double> vec(nmax,0.0); 
  for(G4int i=0; i<nmax; ++i) {
    vec[i] = CrossSection(Z, G4AtomicShellEnumerator(i), kinEnergy); 
  }
  return vec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4double> 
G4LivermoreIonisationCrossSection::Probabilities(G4int Z,
						 G4double kinEnergy,
						 G4double,
						 G4double,
						 const G4Material*)
{
  std::vector<G4double> vec = GetCrossSection(Z, kinEnergy);
  size_t n = vec.size();
  size_t i;
  G4double sum = 0.0;
  for(i=0; i<n; ++i) { sum += vec[i]; }
  if(sum > 0.0) { 
    sum = 1.0/sum; 
    for(i=0; i<n; ++i) { vec[i] = vec[i]*sum; }
  }
  return vec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
