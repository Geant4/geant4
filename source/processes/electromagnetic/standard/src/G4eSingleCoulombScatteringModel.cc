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
//	G4eSingleCoulombScatteringModel.cc
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4eSingleCoulombScatteringModel
//
// Author:      Cristina Consolandi
//
// Creation date: 20.10.2012
//
//	Class Description:
//	Single Scattering model for electron-nuclei interaction.
//	Suitable for high energy electrons and low scattering angles.
//
//
// Reference:
//      M.J. Boschini et al. "Non Ionizing Energy Loss induced by Electrons
//      in the Space Environment" Proc. of the 13th International Conference
//      on Particle Physics and Advanced Technology
//
//	(13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//	Available at: http://arxiv.org/abs/1111.4042v4
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4eSingleCoulombScatteringModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Proton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4UnitsTable.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eSingleCoulombScatteringModel::G4eSingleCoulombScatteringModel(const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(1.0)
{
  fNistManager = G4NistManager::Instance();
  theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  fParticleChange = nullptr;

  pCuts=nullptr;
  currentMaterial = nullptr;
  currentElement  = nullptr;
  currentCouple = nullptr;

  lowEnergyLimit  = 0*keV;
  recoilThreshold = 0.*eV;
  XSectionModel = 1;
  FormFactor = 0;
  particle = nullptr;
  mass=0.0;
  currentMaterialIndex = -1;

  Mottcross = new G4ScreeningMottCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eSingleCoulombScatteringModel::~G4eSingleCoulombScatteringModel()
{
  delete  Mottcross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eSingleCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
						 const G4DataVector&  cuts)
{
  G4EmParameters* param = G4EmParameters::Instance();

  SetupParticle(p);
  currentCouple = nullptr;
  currentMaterialIndex = -1;
  //cosThetaMin = cos(PolarAngleLimit());
  Mottcross->Initialise(p,cosThetaMin);

  pCuts = &cuts;
  //G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);

  /*
  G4cout << "!!! G4eSingleCoulombScatteringModel::Initialise for "
         << part->GetParticleName() << "  cos(TetMin)= " << cosThetaMin
         << "  cos(TetMax)= " << cosThetaMax <<G4endl;
  G4cout << "cut= " << (*pCuts)[0] << "  cut1= " << (*pCuts)[1] << G4endl;
  */

  if(!fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }

  if(IsMaster()) {
    InitialiseElementSelectors(p,cuts);
  }

  FormFactor=param->NuclearFormfactorType();

  //G4cout<<"NUCLEAR FORM FACTOR: "<<FormFactor<<G4endl;
}

void G4eSingleCoulombScatteringModel::InitialiseLocal(const G4ParticleDefinition*,
                                                G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eSingleCoulombScatteringModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy,
				G4double Z,
				G4double ,
				G4double,
				G4double )
{
  SetupParticle(p);

  G4double cross =0.0;
  if(kinEnergy < lowEnergyLimit) return cross;

  DefineMaterial(CurrentCouple());

  //Total Cross section
  Mottcross->SetupKinematic(kinEnergy, Z);
  cross = Mottcross->NuclearCrossSection(FormFactor,XSectionModel);

  //cout<< "Compute Cross Section....cross "<<G4BestUnit(cross,"Surface") << " cm2 "<< cross/cm2 <<" Z: "<<Z<<" kinEnergy: "<<kinEnergy<<endl;

  //G4cout<<"Energy: "<<kinEnergy/MeV<<" Total Cross: "<<cross<<G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eSingleCoulombScatteringModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double cutEnergy,
			       G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  //cout<<"--- kinEnergy "<<kinEnergy<<endl;

  if(kinEnergy < lowEnergyLimit) return;

  DefineMaterial(couple);
  SetupParticle(dp->GetDefinition());

  // Choose nucleus
  //last two :cutEnergy= min e kinEnergy=max
  currentElement = SelectRandomAtom(couple,particle,
				    kinEnergy,cutEnergy,kinEnergy);

  G4double Z  = currentElement->GetZ();
  G4int iz    = G4int(Z);
  G4int ia = SelectIsotopeNumber(currentElement);
  G4double mass2 = G4NucleiProperties::GetNuclearMass(ia, iz);

  //G4cout<<"..Z: "<<Z<<" ..iz: "<<iz<<" ..ia: "<<ia<<" ..mass2: "<<mass2<<G4endl;

  Mottcross->SetupKinematic(kinEnergy, Z);
  G4double cross= Mottcross->NuclearCrossSection(FormFactor,XSectionModel);
  if(cross == 0.0) { return; }
  //cout<< "Energy: "<<kinEnergy/MeV<<" Z: "<<Z<<"....cross "<<G4BestUnit(cross,"Surface") << " cm2 "<< cross/cm2 <<endl;

  G4double z1 = Mottcross->GetScatteringAngle(FormFactor,XSectionModel);
  G4double sint = sin(z1);
  G4double cost = sqrt(1.0 - sint*sint);
  G4double phi  = twopi* G4UniformRand();

  // kinematics in the Lab system
  G4double ptot = dp->GetTotalMomentum();
  G4double e1   = dp->GetTotalEnergy();
  // Lab. system kinematics along projectile direction
  G4LorentzVector v0 = G4LorentzVector(0, 0, ptot, e1);
  G4double bet  = ptot/(v0.e() + mass2);
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));

  //CM Projectile
  G4double momCM = gam*(ptot - bet*e1);
  G4double eCM   = gam*(e1 - bet*ptot);
  //energy & momentum after scattering of incident particle
  G4double pxCM = momCM*sint*cos(phi);
  G4double pyCM = momCM*sint*sin(phi);
  G4double pzCM = momCM*cost;

  //CM--->Lab
  G4LorentzVector v1(pxCM , pyCM, gam*(pzCM + bet*eCM), gam*(eCM + bet*pzCM));

  // Rotate to global system
  G4ThreeVector dir = dp->GetMomentumDirection();
  G4ThreeVector newDirection = v1.vect().unit();
  newDirection.rotateUz(dir);

  fParticleChange->ProposeMomentumDirection(newDirection);

  // recoil
  v0 -= v1;
  G4double trec = v0.e();
  G4double edep = 0.0;

  G4double tcut = recoilThreshold;

  //G4cout<<" Energy Transfered: "<<trec/eV<<G4endl;

  if(pCuts) {
    tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]);
    //G4cout<<"Cuts: "<<(*pCuts)[currentMaterialIndex]/eV<<" eV"<<G4endl;
    //G4cout<<"Threshold: "<<tcut/eV<<" eV"<<G4endl;
  }

  if(trec > tcut) {
    G4ParticleDefinition* ion = theIonTable->GetIon(iz, ia, 0);
    newDirection = v0.vect().unit();
    newDirection.rotateUz(dir);
    G4DynamicParticle* newdp  = new G4DynamicParticle(ion, newDirection, trec);
    fvect->push_back(newdp);
  } else if(trec > 0.0) {
    edep = trec;
    fParticleChange->ProposeNonIonizingEnergyDeposit(edep);
  }

  // finelize primary energy and energy balance
  G4double finalT = v1.e() - mass;
  //G4cout<<"Final Energy: "<<finalT/eV<<G4endl;
  if(finalT <= lowEnergyLimit) {
    edep += finalT;
    finalT = 0.0;
  }
  edep = std::max(edep, 0.0);
  fParticleChange->SetProposedKineticEnergy(finalT);
  fParticleChange->ProposeLocalEnergyDeposit(edep);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
