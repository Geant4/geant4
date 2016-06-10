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
// $Id: G4hCoulombScatteringModel.cc 81579 2014-06-03 10:15:54Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hCoulombScatteringModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 08.06.2012 from G4eCoulombScatteringModel
//
// Modifications:
//
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4hCoulombScatteringModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ElementTable.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Proton.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NucleiProperties.hh"
#include "G4Pow.hh"
#include "G4LossTableManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4hCoulombScatteringModel::G4hCoulombScatteringModel(G4bool combined)
  : G4VEmModel("eCoulombScattering"),
    cosThetaMin(1.0),
    cosThetaMax(-1.0),
    isCombined(combined),
    isInitialised(false)
{
  fParticleChange = 0;
  fNistManager = G4NistManager::Instance();
  theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  theProton   = G4Proton::Proton();
  currentMaterial = 0; 

  pCuts = 0;

  lowEnergyThreshold = 1*keV;  // particle will be killed for lower energy
  recoilThreshold = 0.*keV; // by default does not work

  particle = 0;
  currentCouple = 0;
  wokvi = new G4WentzelVIRelXSection(combined);

  currentMaterialIndex = 0;

  cosTetMinNuc = 1.0;
  cosTetMaxNuc = -1.0;
  elecRatio = 0.0;
  mass = proton_mass_c2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hCoulombScatteringModel::~G4hCoulombScatteringModel()
{
  delete wokvi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector& cuts)
{
  SetupParticle(p);
  currentCouple = 0;

  if(isCombined) {
    cosThetaMin = 1.0;
    G4double tet = PolarAngleLimit();
    if(tet > pi)       { cosThetaMin = -1.0; }
    else if(tet > 0.0) { cosThetaMin = cos(tet); }
  }

  wokvi->Initialise(p, cosThetaMin);
  /*  
  G4cout << "G4hCoulombScatteringModel: " << particle->GetParticleName()
         << "  1-cos(ThetaLimit)= " << 1 - cosThetaMin
	 << "  cos(thetaMax)= " <<  cosThetaMax
	 << G4endl;
  */
  pCuts = &cuts;
  //  G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);
  //G4cout << "!!! G4hCoulombScatteringModel::Initialise for " 
  //	 << p->GetParticleName() << "  cos(TetMin)= " << cosThetaMin 
  //	 << "  cos(TetMax)= " << cosThetaMax <<G4endl;
  // G4cout << "cut0= " << cuts[0] << "  cut1= " << cuts[1] << G4endl;
  if(!isInitialised) {
    isInitialised = true;
    fParticleChange = GetParticleChangeForGamma();
  }
  if(mass < GeV) {
    InitialiseElementSelectors(p,cuts);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hCoulombScatteringModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double,
		G4double cutEnergy, G4double)
{
  //G4cout << "### G4hCoulombScatteringModel::ComputeCrossSectionPerAtom  for " 
  //  << p->GetParticleName()<<" Z= "<<Z<<" e(MeV)= "<< kinEnergy/MeV << G4endl; 
  G4double cross = 0.0;
  if(p != particle) { SetupParticle(p); }

  // cross section is set to zero to avoid problems in sample secondary
  if(kinEnergy <= 0.0) { return cross; }
  DefineMaterial(CurrentCouple());
  cosTetMinNuc = wokvi->SetupKinematic(kinEnergy, currentMaterial);
  if(cosThetaMax < cosTetMinNuc) {
    G4int iz = G4int(Z);
    cosTetMinNuc = wokvi->SetupTarget(iz, cutEnergy);
    cosTetMaxNuc = cosThetaMax; 
    if(iz == 1 && cosTetMaxNuc < 0.0 && particle == theProton) { 
      cosTetMaxNuc = 0.0; 
    }
    cross =  wokvi->ComputeNuclearCrossSection(cosTetMinNuc, cosTetMaxNuc);
    elecRatio = wokvi->ComputeElectronCrossSection(cosTetMinNuc, cosThetaMax);
    cross += elecRatio;
    if(cross > 0.0) { elecRatio /= cross; }  
  }
  /*
  if(p->GetParticleName() == "e-") 
  G4cout << "e(MeV)= " << kinEnergy/MeV << " cross(b)= " << cross/barn  
	 << " 1-cosTetMinNuc= " << 1-cosTetMinNuc
	 << " 1-cosTetMaxNuc= " << 1-cosTetMaxNuc
	 << G4endl;
  */
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hCoulombScatteringModel::SampleSecondaries(
                std::vector<G4DynamicParticle*>* fvect,
		const G4MaterialCutsCouple* couple,
		const G4DynamicParticle* dp,
		G4double cutEnergy,
		G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();

  // absorb particle below low-energy limit to avoid situation
  // when a particle has no energy loss
  if(kinEnergy < lowEnergyThreshold) { 
    fParticleChange->SetProposedKineticEnergy(0.0);
    fParticleChange->ProposeLocalEnergyDeposit(kinEnergy);
    fParticleChange->ProposeNonIonizingEnergyDeposit(kinEnergy);
    return; 
  }
  SetupParticle(dp->GetDefinition());
  DefineMaterial(couple);

  //G4cout << "G4hCoulombScatteringModel::SampleSecondaries e(MeV)= " 
  //	 << kinEnergy << "  " << particle->GetParticleName() 
  //	 << " cut= " << cutEnergy<< G4endl;
 
  // Choose nucleus
  const G4Element* currentElement = 
    SelectRandomAtom(couple,particle,kinEnergy,cutEnergy,kinEnergy);

  G4double Z = currentElement->GetZ();

  if(ComputeCrossSectionPerAtom(particle,kinEnergy, Z,
  				kinEnergy, cutEnergy, kinEnergy) == 0.0) 
    { return; }

  G4int iz = G4int(Z);
  G4int ia = SelectIsotopeNumber(currentElement);
  G4double targetMass = G4NucleiProperties::GetNuclearMass(ia, iz);
  wokvi->SetTargetMass(targetMass);

  G4ThreeVector newDirection = 
    wokvi->SampleSingleScattering(cosTetMinNuc, cosThetaMax, elecRatio);
  G4double cost = newDirection.z();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   

  // recoil sampling assuming a small recoil
  // and first order correction to primary 4-momentum
  G4double mom2 = wokvi->GetMomentumSquare();
  G4double trec = mom2*(1.0 - cost)
    /(targetMass + (mass + kinEnergy)*(1.0 - cost));

  // the check likely not needed
  if(trec > kinEnergy) { trec = kinEnergy; }
  G4double finalT = kinEnergy - trec; 
  G4double edep = 0.0;
  //G4cout<<"G4eCoulombScatteringModel: finalT= "<<finalT<<" Trec= "
  //	<<trec << " Z= " << iz << " A= " << ia<<G4endl;

  G4double tcut = recoilThreshold;
  if(pCuts) { tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]); }

  if(trec > tcut) {
    G4ParticleDefinition* ion = theIonTable->GetIon(iz, ia, 0);
    G4ThreeVector dir = (direction*sqrt(mom2) - 
			 newDirection*sqrt(finalT*(2*mass + finalT))).unit();
    G4DynamicParticle* newdp = new G4DynamicParticle(ion, dir, trec);
    fvect->push_back(newdp);
  } else {
    edep = trec;
    fParticleChange->ProposeNonIonizingEnergyDeposit(edep);
  }

  // finelize primary energy and energy balance
  // this threshold may be applied only because for low-enegry
  // e+e- msc model is applied
  if(finalT <= lowEnergyThreshold) { 
    edep += finalT;  
    finalT = 0.0;
  } 
  fParticleChange->SetProposedKineticEnergy(finalT);
  fParticleChange->ProposeLocalEnergyDeposit(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


