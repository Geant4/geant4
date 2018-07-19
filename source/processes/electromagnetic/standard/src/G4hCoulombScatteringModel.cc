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
// $Id: G4hCoulombScatteringModel.cc 104802 2017-06-19 07:11:40Z gcosmo $
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
  : G4VEmModel("hCoulombScattering"),
    cosThetaMin(1.0),
    cosThetaMax(-1.0),
    isCombined(combined)
{
  fParticleChange = nullptr;
  fNistManager = G4NistManager::Instance();
  theIonTable  = G4ParticleTable::GetParticleTable()->GetIonTable();
  theProton    = G4Proton::Proton();
  currentMaterial = nullptr; 
  fixedCut = -1.0;

  pCuts = nullptr;

  recoilThreshold = 0.0; // by default does not work

  particle = nullptr;
  currentCouple = nullptr;
  wokvi = new G4WentzelVIRelXSection();

  currentMaterialIndex = 0;
  mass = CLHEP::proton_mass_c2;
  elecRatio = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4hCoulombScatteringModel::~G4hCoulombScatteringModel()
{
  delete wokvi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4hCoulombScatteringModel::Initialise(const G4ParticleDefinition* part,
					   const G4DataVector& cuts)
{
  SetupParticle(part);
  currentCouple = nullptr;

  // defined theta limit between single and multiple scattering 
  isCombined = true;
  G4double tet = PolarAngleLimit();

  if(tet <= 0.0) { 
    cosThetaMin = 1.0; 
    isCombined = false; 
  } else if(tet >= CLHEP::pi) { 
    cosThetaMin = -1.0; 
  } else { 
    cosThetaMin = cos(tet); 
  }

  wokvi->Initialise(part, cosThetaMin);
  /*
  G4cout << "G4hCoulombScatteringModel: " << particle->GetParticleName()
         << "  1-cos(ThetaLimit)= " << 1 - cosThetaMin
	 << "  cos(thetaMax)= " <<  cosThetaMax
	 << G4endl;
  */
  pCuts = &cuts;
  //G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);
  /*
  G4cout << "!!! G4hCoulombScatteringModel::Initialise for " 
  	 << part->GetParticleName() << "  cos(TetMin)= " << cosThetaMin 
  	 << "  cos(TetMax)= " << cosThetaMax <<G4endl;
  G4cout << "cut= " << (*pCuts)[0] << "  cut1= " << (*pCuts)[1] << G4endl;
  */
  if(!fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
  }
  if(IsMaster() && mass < GeV && part->GetParticleName() != "GenericIon") {
    InitialiseElementSelectors(part, cuts);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4hCoulombScatteringModel::InitialiseLocal(const G4ParticleDefinition*, 
						G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4hCoulombScatteringModel::MinPrimaryEnergy(const G4Material* material,
					    const G4ParticleDefinition* part,
					    G4double)
{
  SetupParticle(part);

  // define cut using cuts for proton
  G4double cut = 
    std::max(recoilThreshold, (*pCuts)[CurrentCouple()->GetIndex()]);

  // find out lightest element
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nelm = material->GetNumberOfElements();

  // select lightest element
  G4int Z = 300;
  for (G4int j=0; j<nelm; ++j) {
    Z = std::min(Z,(*theElementVector)[j]->GetZasInt());
  }
  G4int A = G4lrint(fNistManager->GetAtomicMassAmu(Z));
  G4double targetMass = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double t = std::max(cut, 0.5*(cut + sqrt(2*cut*targetMass)));

  return t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4hCoulombScatteringModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double,
		G4double cutEnergy, G4double)
{
  //G4cout << "### G4hCoulombScatteringModel::ComputeCrossSectionPerAtom  for " 
  //<< p->GetParticleName()<<" Z= "<<Z<<" e(MeV)= "<< kinEnergy/MeV << G4endl; 
  G4double cross = 0.0;
  elecRatio = 0.0;
  if(p != particle) { SetupParticle(p); }

  // cross section is set to zero to avoid problems in sample secondary
  if(kinEnergy <= 0.0) { return cross; }
  DefineMaterial(CurrentCouple());

  G4int iz = G4lrint(Z);
  G4double tmass = (1 == iz) ? proton_mass_c2 :
    fNistManager->GetAtomicMassAmu(iz)*amu_c2; 
  wokvi->SetTargetMass(tmass);

  G4double costmin = 
    wokvi->SetupKinematic(kinEnergy, currentMaterial);

  if(cosThetaMax < costmin) {
    G4double cut = (0.0 < fixedCut) ? fixedCut : cutEnergy;
    costmin = wokvi->SetupTarget(iz, cut);
    G4double costmax = 
      (1 == iz && particle == theProton && cosThetaMax < 0.0) 
      ? 0.0 : cosThetaMax; 
    if(costmin > costmax) {
      cross = wokvi->ComputeNuclearCrossSection(costmin, costmax)
	+ wokvi->ComputeElectronCrossSection(costmin, costmax);
    }
    /*  
  if(p->GetParticleName() == "mu+") 
  G4cout << "e(MeV)= " << kinEnergy/MeV << " cross(b)= " << cross/barn  
	 << " 1-costmin= " << 1-costmin
	 << " 1-costmax= " << 1-costmax
	 << " 1-cosThetaMax= " << 1-cosThetaMax
	 << G4endl;
    */
  }
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
  SetupParticle(dp->GetDefinition());
  DefineMaterial(couple);

  // Choose nucleus
  G4double cut = (0.0 < fixedCut) ? fixedCut : cutEnergy;

  const G4Element* elm = SelectRandomAtom(couple,particle,
					  kinEnergy,cut,kinEnergy);

  G4int iz = elm->GetZasInt();
  G4int ia = SelectIsotopeNumber(elm);
  G4double mass2 = G4NucleiProperties::GetNuclearMass(ia, iz);

  wokvi->SetTargetMass(mass2);
  wokvi->SetupKinematic(kinEnergy, currentMaterial);
  G4double costmin = wokvi->SetupTarget(iz, cut);
  G4double costmax = (1 == iz && particle == theProton && cosThetaMax < 0.0) 
    ? 0.0 :  cosThetaMax; 
  if(costmin <= costmax) { return; }

  G4double cross = wokvi->ComputeNuclearCrossSection(costmin, costmax);
  G4double ecross = wokvi->ComputeElectronCrossSection(costmin, costmax);
  G4double ratio = ecross/(cross + ecross);

  G4ThreeVector newDirection = 
    wokvi->SampleSingleScattering(costmin, costmax, ratio);

  // kinematics in the Lab system
  G4double ptot = dp->GetTotalMomentum();
  G4double e1   = dp->GetTotalEnergy();
  
  // Lab. system kinematics along projectile direction
  G4LorentzVector v0 = G4LorentzVector(0, 0, ptot, e1 + mass2);
  G4double bet  = ptot/v0.e();
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));

  // CM projectile
  G4double momCM = gam*(ptot - bet*e1); 
  G4double eCM   = gam*(e1 - bet*ptot); 
  // energy & momentum after scattering of incident particle
  G4double pxCM = momCM*newDirection.x();
  G4double pyCM = momCM*newDirection.y();
  G4double pzCM = momCM*newDirection.z();

  // CM--->Lab
  G4LorentzVector v1(pxCM , pyCM, gam*(pzCM + bet*eCM), gam*(eCM + bet*pzCM));

  G4ThreeVector dir = dp->GetMomentumDirection(); 
  newDirection = v1.vect().unit();
  newDirection.rotateUz(dir);   
  
  fParticleChange->ProposeMomentumDirection(newDirection);   
  
  // recoil
  v0 -= v1; 
  G4double trec = v0.e() - mass2;
  G4double edep = 0.0;

  G4double tcut = recoilThreshold;
  if(pCuts) { tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]); }
 
  if(trec > tcut) {
    G4ParticleDefinition* ion = theIonTable->GetIon(iz, ia, 0);
    newDirection = v0.vect().unit();
    newDirection.rotateUz(dir);   
    G4DynamicParticle* newdp = new G4DynamicParticle(ion, newDirection, trec);
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
  fParticleChange->SetProposedKineticEnergy(finalT);
  fParticleChange->ProposeLocalEnergyDeposit(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
