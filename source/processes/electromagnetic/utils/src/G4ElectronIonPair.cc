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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ElectronIonPair
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 08.07.2008
//
// Modifications:
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ElectronIonPair.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4Track.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ElectronIonPair::G4ElectronIonPair(G4int verb)
{
  verbose = verb;
  curMaterial = nullptr;
  curMeanEnergy = 0.0;
  nMaterials = 0;
  invFanoFactor = 1.0/0.2;
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ElectronIonPair::~G4ElectronIonPair() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ElectronIonPair::MeanNumberOfIonsAlongStep(
				   const G4ParticleDefinition* part, 
				   const G4Material* material,
				   G4double edep,
				   G4double niel)
{
  G4double nion = 0.0;

  // NIEL does not provide ionisation clusters
  if(edep > niel) {

    // neutral particles do not produce ionisation along step
    if(part->GetPDGCharge() != 0.0) {

      // select material
      if(material != curMaterial) {
	curMaterial = material;
	curMeanEnergy = material->GetIonisation()->GetMeanEnergyPerIonPair();

	// if mean energy is not defined then look into G4 DB
	if(0.0 == curMeanEnergy) {
	  curMeanEnergy = FindG4MeanEnergyPerIonPair(material);
	} 
      }
      if(curMeanEnergy > 0.0) { nion = (edep - niel)/curMeanEnergy; }
    }
  }
  return nion;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4ThreeVector>* 
G4ElectronIonPair::SampleIonsAlongStep(const G4Step* step)
{
  std::vector<G4ThreeVector>* v = nullptr;

  G4int nion = SampleNumberOfIonsAlongStep(step);

  // sample ionisation along step
  if(nion > 0) {

    v = new std::vector<G4ThreeVector>;
    G4ThreeVector prePos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector deltaPos = step->GetPostStepPoint()->GetPosition() - prePos;  
    for(G4int i=0; i<nion; ++i) {
      v->push_back( prePos + deltaPos*G4UniformRand() );
    }
    if(verbose > 1 ) { 
      G4cout << "### G4ElectronIonPair::SampleIonisationPoints: "
	     << v->size() << "  ion pairs are added" << G4endl;
    }
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int 
G4ElectronIonPair::ResidualeChargePostStep(const G4ParticleDefinition*,
					   const G4TrackVector*,
					   G4int subType) const
{
  G4int nholes = 0;

  if(2 == subType || 12 == subType || 13 == subType) { nholes = 1; }
  return nholes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ElectronIonPair::FindG4MeanEnergyPerIonPair(const G4Material* mat) const
{
  G4String name = mat->GetName();
  G4double res  = 0.0;

  // is this material in the vector?  
  for(G4int j=0; j<nMaterials; ++j) {
    if(name == g4MatNames[j]) {
      res = g4MatData[j];
      mat->GetIonisation()->SetMeanEnergyPerIonPair(res);
      if(verbose > 0) {
	G4cout << "### G4ElectronIonPair::FindG4MeanEnergyPerIonPair for "
	       << name << " Epair= " << res/eV << " eV is set"
	       << G4endl;
      }
      break;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ElectronIonPair:: DumpMeanEnergyPerIonPair() const
{
  std::size_t nmat = G4Material::GetNumberOfMaterials();
  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  if(nmat > 0) {
    G4cout << "### G4ElectronIonPair: mean energy per ion pair available:" 
	   << G4endl;
    for(std::size_t i=0; i<nmat; ++i) {
      const G4Material* mat = (*mtable)[i];
      G4double x = mat->GetIonisation()->GetMeanEnergyPerIonPair();
      if(x > 0.0) {
	G4cout << "   " << mat->GetName() << "   Epair=  " 
	       << x/eV << " eV" << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ElectronIonPair::DumpG4MeanEnergyPerIonPair() const
{
  if(nMaterials > 0) {
    G4cout << "### G4ElectronIonPair: mean energy per ion pair "
	   << " for Geant4 materials" << G4endl;
    for(G4int i=0; i<nMaterials; ++i) {
      G4cout << "   " << g4MatNames[i] << "    Epair= " 
	     << g4MatData[i]/eV << " eV" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ElectronIonPair::Initialise()
{
  // ICRU Report 31, 1979
  g4MatNames.push_back("G4_Si");
  g4MatData.push_back(3.62*eV);

  g4MatNames.push_back("G4_Ge");
  g4MatData.push_back(2.97*eV);

  g4MatNames.push_back("G4_He");
  g4MatData.push_back(44.4*eV);

  g4MatNames.push_back("G4_N");
  g4MatData.push_back(36.4*eV);

  g4MatNames.push_back("G4_O");
  g4MatData.push_back(32.3*eV);

  g4MatNames.push_back("G4_Ne");
  g4MatData.push_back(36.8*eV);

  g4MatNames.push_back("G4_Ar");
  g4MatData.push_back(26.34*eV);

  g4MatNames.push_back("G4_Kr");
  g4MatData.push_back(24.1*eV);

  g4MatNames.push_back("G4_Xe");
  g4MatData.push_back(21.6*eV);

  g4MatNames.push_back("G4_lAr");
  g4MatData.push_back(23.6*eV);

  g4MatNames.push_back("G4_lKr");
  g4MatData.push_back(20.5*eV);

  g4MatNames.push_back("G4_lXe");
  g4MatData.push_back(15.6*eV);

  g4MatNames.push_back("G4_AIR");
  g4MatData.push_back(35.1*eV);

  nMaterials = (G4int)g4MatData.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
