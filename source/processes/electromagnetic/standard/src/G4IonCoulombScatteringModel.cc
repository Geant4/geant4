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
//	G4IonCoulombScatteringModel.cc
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4IonCoulombScatteringModel
//
// Author:      Cristina Consolandi
//
// Creation date: 05.10.2010 from G4eCoulombScatteringModel 
//                               & G4CoulombScatteringModel
//
// Class Description:
//      Single Scattering Model for
//      for protons, alpha and heavy Ions
//
// Reference:
//      M.J. Boschini et al. "Nuclear and Non-Ionizing Energy-Loss 
//      for Coulomb ScatteredParticles from Low Energy up to Relativistic 
//      Regime in Space Radiation Environment"
//      Accepted for publication in the Proceedings of  the  ICATPP Conference
//      on Cosmic Rays for Particle and Astroparticle Physics, Villa  Olmo, 7-8
//      October,  2010, to be published by World Scientific (Singapore).
//
//      Available for downloading at:
//	http://arxiv.org/abs/1011.4822
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4IonCoulombScatteringModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Proton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4IonCoulombScatteringModel::G4IonCoulombScatteringModel(const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(1.0)
{
  fNistManager = G4NistManager::Instance();
  theIonTable  = G4ParticleTable::GetParticleTable()->GetIonTable();
  theProton    = G4Proton::Proton();

  pCuts = nullptr;
  currentMaterial = nullptr;
  currentElement  = nullptr;
  currentCouple   = nullptr;
  fParticleChange = nullptr;

  recoilThreshold = 0.*eV;
  heavycorr =0;
  particle = nullptr;
  mass=0;
  currentMaterialIndex = -1;

  ioncross = new G4IonCoulombCrossSection(); 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonCoulombScatteringModel::~G4IonCoulombScatteringModel()
{ 
  delete  ioncross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
					     const G4DataVector& cuts)
{
  SetupParticle(p);
  currentCouple = nullptr;
  currentMaterialIndex = -1;
  ioncross->Initialise(p,cosThetaMin);
 
  pCuts = &cuts;
  //  G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);
  if(!fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4IonCoulombScatteringModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy, 
				G4double Z, 
				G4double, G4double, G4double)
{
  SetupParticle(p);
 
  G4double cross = 0.0;

  DefineMaterial(CurrentCouple());

  G4int iz = G4lrint(Z);

  //from lab to pCM & mu_rel of effective particle
  G4double tmass = proton_mass_c2;
  if(1 < iz) {
    tmass = fNistManager->GetAtomicMassAmu(iz)*amu_c2;
  }
  ioncross->SetupKinematic(kinEnergy, tmass);
  ioncross->SetupTarget(Z, kinEnergy, heavycorr);
  cross = ioncross->NuclearCrossSection();
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombScatteringModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double, G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  DefineMaterial(couple);
  SetupParticle(dp->GetDefinition());

  // Choose nucleus
  currentElement = SelectTargetAtom(couple, particle, kinEnergy,
                                    dp->GetLogKineticEnergy());

  G4int iz = currentElement->GetZasInt();
  G4int ia = SelectIsotopeNumber(currentElement);
  G4double mass2 = G4NucleiProperties::GetNuclearMass(ia, iz);

  ioncross->SetupKinematic(kinEnergy, mass2);
  ioncross->SetupTarget(currentElement->GetZ(), kinEnergy, heavycorr);
    
  //scattering angle, z1 == (1-cost)
  G4double z1 = ioncross->SampleCosineTheta(); 
  if(z1 > 2.0)      { z1 = 2.0; }
  else if(z1 < 0.0) { z1 = 0.0; }
  /*
  G4cout << "Sample: " << particle->GetParticleName() 
	 << " mass(GeV)= " << mass/GeV 
	 << " Ekin(MeV)= " << kinEnergy << " cost= " << 1. - z1 << G4endl; 
  G4cout << "     Z= " << iz << " A= " << ia 
	 << " mass(GeV)= " << mass2/GeV << G4endl;
  */
  G4double cost = 1.0 - z1;
  G4double sint = sqrt(z1*(1.0 + cost));
  G4double phi  = twopi * G4UniformRand();

  // kinematics in the Lab system
  G4double ptot = sqrt(kinEnergy*(kinEnergy + 2.0*mass));
  G4double e1   = mass + kinEnergy;
  
  // Lab. system kinematics along projectile direction
  G4LorentzVector v0 = G4LorentzVector(0, 0, ptot, e1+mass2);
  G4LorentzVector v1 = G4LorentzVector(0, 0, ptot, e1);
  G4ThreeVector bst = v0.boostVector();
  v1.boost(-bst);
  // CM projectile
  G4double momCM = v1.pz(); 
  
  // Momentum after scattering of incident particle
  v1.setX(momCM*sint*cos(phi));
  v1.setY(momCM*sint*sin(phi));
  v1.setZ(momCM*cost);

  // CM--->Lab
  v1.boost(bst);

  // Rotate to global system
  G4ThreeVector dir = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection = v1.vect().unit();
  newDirection.rotateUz(dir);   
  
  fParticleChange->ProposeMomentumDirection(newDirection);   
  
  // recoil v0 energy is kinetic
  v0 -= v1; 
  G4double trec = std::max(v0.e() - mass2, 0.0);
  G4double edep = 0.0;

  G4double tcut = recoilThreshold;
  if(pCuts) { 
    tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]); 
    //G4cout<<" tcut eV "<<tcut/eV<<endl;
  }
 
  // Recoil
  if(trec > tcut) {
    G4ParticleDefinition* ion = theIonTable->GetIon(iz, ia, 0);
    newDirection = v0.vect().unit();
    newDirection.rotateUz(dir);   
    auto newdp = new G4DynamicParticle(ion, newDirection, trec);
    fvect->push_back(newdp);
  } else if(trec > 0.0) {
    edep = trec;
    fParticleChange->ProposeNonIonizingEnergyDeposit(edep);
  }

  // finelize primary energy and energy balance
  G4double finalT = v1.e() - mass;
  if(finalT < 0.0) { 
    edep += finalT;
    finalT = 0.0;
  } 
  edep = std::max(edep, 0.0);
  //G4cout << "Efinal(MeV)= " << finalT << " Edep(MeV)= " << edep 
  //	 << " Trec(MeV)= " << trec << G4endl;
  fParticleChange->SetProposedKineticEnergy(finalT);
  fParticleChange->ProposeLocalEnergyDeposit(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
		
