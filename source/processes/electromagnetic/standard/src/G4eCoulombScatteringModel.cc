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
// $Id: G4eCoulombScatteringModel.cc,v 1.91 2010-11-13 18:45:55 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eCoulombScatteringModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 22.08.2005
//
// Modifications:
//
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.08.06 V.Ivanchenko add inline function ScreeningParameter 
// 09.10.07 V.Ivanchenko reorganized methods, add cut dependence in scattering off e- 
// 09.06.08 V.Ivanchenko add SelectIsotope and sampling of the recoil ion 
// 16.06.09 C.Consolandi fixed computation of effective mass
// 27.05.10 V.Ivanchenko added G4WentzelOKandVIxSection class to
//              compute cross sections and sample scattering angle
//
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eCoulombScatteringModel.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ElementTable.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Proton.hh"
#include "G4ParticleTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NucleiProperties.hh"
#include "G4Pow.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eCoulombScatteringModel::G4eCoulombScatteringModel(const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(1.0),
    cosThetaMax(-1.0),
    isInitialised(false)
{
  fParticleChange = 0; 
  fNistManager = G4NistManager::Instance();
  theParticleTable = G4ParticleTable::GetParticleTable();
  theProton   = G4Proton::Proton();
  currentMaterial = 0; 
  currentElement  = 0;

  pCuts = 0;

  lowEnergyLimit = 1*keV;
  recoilThreshold = 0.*keV;
  particle = 0;
  currentCouple = 0;
  wokvi = new G4WentzelOKandVIxSection();

  currentMaterialIndex = 0;

  cosTetMinNuc = 1.0;
  cosTetMaxNuc = -1.0;
  elecRatio = 0.0;
  mass = proton_mass_c2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eCoulombScatteringModel::~G4eCoulombScatteringModel()
{
  delete wokvi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector& cuts)
{
  SetupParticle(p);
  currentCouple = 0;
  cosThetaMin = cos(PolarAngleLimit());
  wokvi->Initialise(p, cosThetaMin);
  /*
  G4cout << "G4eCoulombScatteringModel: factorA2(GeV^2) = " << factorA2/(GeV*GeV) 
         << "  1-cos(ThetaLimit)= " << 1 - cosThetaMin
	 << "  cos(thetaMax)= " <<  cosThetaMax
	 << G4endl;
  */
  pCuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);
  //G4cout << "!!! G4eCoulombScatteringModel::Initialise for " 
  //	 << p->GetParticleName() << "  cos(TetMin)= " << cosThetaMin 
  //	 << "  cos(TetMax)= " << cosThetaMax <<G4endl;
  // G4cout << "cut0= " << cuts[0] << "  cut1= " << cuts[1] << G4endl;
  if(!isInitialised) {
    isInitialised = true;
    fParticleChange = GetParticleChangeForGamma();
  }
  if(mass < GeV && particle->GetParticleType() != "nucleus") {
    InitialiseElementSelectors(p,cuts);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double,
		G4double cutEnergy, G4double)
{
  //G4cout << "### G4eCoulombScatteringModel::ComputeCrossSectionPerAtom  for " 
  //  << p->GetParticleName()<<" Z= "<<Z<<" e(MeV)= "<< kinEnergy/MeV << G4endl; 
  G4double xsec = 0.0;
  if(p != particle) { SetupParticle(p); }

  // cross section is set to zero to avoid problems in sample secondary
  if(kinEnergy <= 0.0) { return xsec; }
  DefineMaterial(CurrentCouple());
  cosTetMinNuc = wokvi->SetupKinematic(kinEnergy, currentMaterial);
  if(cosThetaMax < cosTetMinNuc) {
    G4int iz = G4int(Z);
    cosTetMinNuc = wokvi->SetupTarget(iz, cutEnergy);
    cosTetMaxNuc = cosThetaMax; 
    if(iz == 1 && cosTetMaxNuc < 0.0 && particle == theProton) { 
      cosTetMaxNuc = 0.0; 
    }
    xsec =  wokvi->ComputeNuclearCrossSection(cosTetMinNuc, cosTetMaxNuc);
    elecRatio = wokvi->ComputeElectronCrossSection(cosTetMinNuc, cosThetaMax);
    xsec += elecRatio;
    if(xsec > 0.0) { elecRatio /= xsec; }  
  }
  /*
  G4cout << "e(MeV)= " << kinEnergy/MeV << " xsec(b)= " << xsec/barn  
	 << " 1-cosTetMinNuc= " << 1-cosTetMinNuc
	 << " 1-cosTetMaxNuc2= " << 1-cosTetMaxNuc2
	 << " 1-cosTetMaxElec= " << 1-cosTetMaxElec
	 << " screenZ= " << screenZ
	 << " formfactA= " << formfactA << G4endl;
  */
  return xsec;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::SampleSecondaries(
                std::vector<G4DynamicParticle*>* fvect,
		const G4MaterialCutsCouple* couple,
		const G4DynamicParticle* dp,
		G4double cutEnergy,
		G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  if(kinEnergy < lowEnergyLimit) {
    fParticleChange->SetProposedKineticEnergy(0.0);
    fParticleChange->ProposeLocalEnergyDeposit(kinEnergy);
    fParticleChange->ProposeNonIonizingEnergyDeposit(kinEnergy);
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange->ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange->ProposeTrackStatus(fStopAndKill); }
    return;
  }
  SetupParticle(dp->GetDefinition());
  DefineMaterial(couple);

  //G4cout << "G4eCoulombScatteringModel::SampleSecondaries e(MeV)= " 
  //	 << kinEnergy << "  " << particle->GetParticleName() 
  //	 << " cut= " << cutEnergy<< G4endl;
 
  // Choose nucleus
  currentElement = SelectRandomAtom(couple,particle,
				    kinEnergy,cutEnergy,kinEnergy);

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
  G4double trec = mom2*(1.0 - cost)/(targetMass + (mass + kinEnergy)*(1.0 - cost));
  G4double finalT = kinEnergy - trec; 
  //G4cout<<"G4eCoulombScatteringModel: finalT= "<<finalT<<" Trec= "<<trec<<G4endl;
  if(finalT <= lowEnergyLimit) { 
    trec = kinEnergy;  
    finalT = 0.0;
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange->ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange->ProposeTrackStatus(fStopAndKill); }
  } 

  fParticleChange->SetProposedKineticEnergy(finalT);
  G4double tcut = recoilThreshold;
  if(pCuts) { tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]); }

  if(trec > tcut) {
    G4ParticleDefinition* ion = theParticleTable->GetIon(iz, ia, 0.0);
    G4ThreeVector dir = (direction*sqrt(mom2) - 
			 newDirection*sqrt(finalT*(2*mass + finalT))).unit();
    G4DynamicParticle* newdp = new G4DynamicParticle(ion, dir, trec);
    fvect->push_back(newdp);
  } else {
    fParticleChange->ProposeLocalEnergyDeposit(trec);
    fParticleChange->ProposeNonIonizingEnergyDeposit(trec);
  }
 
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


