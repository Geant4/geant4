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
// $Id: G4eCoulombScatteringModel.cc,v 1.59 2008/10/22 18:39:29 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.08.06 V.Ivanchenko add inline function ScreeningParameter 
// 09.10.07 V.Ivanchenko reorganized methods, add cut dependence in scattering off e- 
// 09.06.08 V.Ivanchenko add SelectIsotope and sampling of the recoil ion 
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
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eCoulombScatteringModel::G4eCoulombScatteringModel(const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(1.0),
    cosThetaMax(-1.0),
    q2Limit(TeV*TeV),
    alpha2(fine_structure_const*fine_structure_const),
    faclim(100.0),
    isInitialised(false)
{
  fNistManager = G4NistManager::Instance();
  theParticleTable = G4ParticleTable::GetParticleTable();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  theProton   = G4Proton::Proton();
  currentMaterial = 0; 
  currentElement  = 0;
  a0 = alpha2*electron_mass_c2*electron_mass_c2/(0.885*0.885);
  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff  = twopi*p0*p0;
  constn = 6.937e-6/(MeV*MeV);
  tkin = targetZ = mom2 = DBL_MIN;
  elecXSection = nucXSection = 0.0;
  recoilThreshold = DBL_MAX;
  ecut = DBL_MAX;
  particle = 0;
  currentCouple = 0;
  for(size_t j=0; j<100; j++) {
    FF[j] = 0.0;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eCoulombScatteringModel::~G4eCoulombScatteringModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector& cuts)
{
  SetupParticle(p);
  currentCouple = 0;
  elecXSection = nucXSection = 0.0;
  tkin = targetZ = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  cosThetaMin = cos(PolarAngleLimit());
  currentCuts = &cuts;
  //G4cout << "!!! G4eCoulombScatteringModel::Initialise for " 
  //	 << p->GetParticleName() << "  cos(TetMin)= " << cosThetaMin 
  //	 << "  cos(TetMax)= " << cosThetaMax <<G4endl;
  if(!isInitialised) {
    isInitialised = true;

    if(pParticleChange)
      fParticleChange = 
	reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChange = new G4ParticleChangeForGamma();
  }
  if(mass < GeV && particle->GetParticleType() != "nucleus") {
    InitialiseElementSelectors(p,cuts);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4eCoulombScatteringModel::ComputeMaxElectronScattering(G4double cutEnergy)
{
  ecut = cutEnergy;
  G4double tmax = tkin;
  cosTetMaxElec = 1.0;
  if(mass > MeV) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = tkin/mass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
    cosTetMaxElec = 1.0 - std::min(cutEnergy, tmax)*electron_mass_c2/mom2;
  } else {

    if(particle == theElectron) tmax *= 0.5;
    G4double t = std::min(cutEnergy, tmax);
    G4double mom21 = t*(t + 2.0*electron_mass_c2);
    G4double t1 = tkin - t;
    //G4cout << "tkin= " << tkin << " t= " << t << " t1= " << t1 << G4endl;
    if(t1 > 0.0) {
      G4double mom22 = t1*(t1 + 2.0*mass);
      G4double ctm = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
      //G4cout << "ctm= " << ctm << G4endl;
      if(ctm < 1.0) cosTetMaxElec = ctm;
    }
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
  SetupParticle(p);
  G4double ekin = std::max(lowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  if(cosTetMaxNuc < cosTetMinNuc) {
    SetupTarget(Z, ekin);
    xsec = CrossSectionPerAtom();  
  }
  /*
  G4cout << "e(MeV)= " << ekin/MeV << "cosTetMinNuc= " << cosTetMinNuc
	 << " cosTetMaxNuc= " << cosTetMaxNuc
	 << " cosTetMaxElec= " << cosTetMaxElec
	 << " screenZ= " << screenZ
	 << " formfactA= " << formfactA
	 << " cosTetMaxHad= " << cosTetMaxHad << G4endl;
  */
  return xsec;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::CrossSectionPerAtom()
{
  // This method needs initialisation before be called

  G4double fac = coeff*targetZ*chargeSquare*invbeta2/mom2; 
  elecXSection = 0.0;
  nucXSection  = 0.0;

  G4double x  = 1.0 - cosTetMinNuc;
  G4double x1 = x + screenZ;

  if(cosTetMaxElec2 < cosTetMinNuc) {
    elecXSection = fac*(cosTetMinNuc - cosTetMaxElec2)/
      (x1*(1.0 - cosTetMaxElec2 + screenZ));
    nucXSection  = elecXSection;
  }

  //G4cout << "XS tkin(MeV)= " << tkin<<" xs= " <<nucXSection 
  //	 << " costmax= " << cosTetMaxNuc2 
  //	 << " costmin= " << cosTetMinNuc << "  Z= " << targetZ <<G4endl;
  if(cosTetMaxNuc2 < cosTetMinNuc) {
    G4double s  = screenZ*formfactA;
    G4double z1 = 1.0 - cosTetMaxNuc2 + screenZ;
    G4double d  = (1.0 - s)/formfactA;
    //G4cout <<"x1= "<<x1<<" z1= " <<z1<<" s= "<<s << " d= " <<d <<G4endl;
    if(d < 0.2*x1) {
      G4double x2 = x1*x1;
      G4double z2 = z1*z1;
      x = (1.0/(x1*x2) - 1.0/(z1*z2) - d*1.5*(1.0/(x2*x2) - 1.0/(z2*z2)))/
	(3.0*formfactA*formfactA);
    } else {
      G4double x2 = x1 + d;
      G4double z2 = z1 + d;
      x = (1.0 + 2.0*s)*((cosTetMinNuc - cosTetMaxNuc2)*(1.0/(x1*z1) + 1.0/(x2*z2)) -
			 2.0*log(z1*x2/(z2*x1))/d);
    }
    nucXSection += fac*targetZ*x;
  }

  //G4cout<<" cross(bn)= "<<nucXSection/barn<<" xsElec(bn)= "<<elecXSection/barn
  //	<< " Asc= " << screenZ << G4endl; 
  
  return nucXSection;
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
  if(kinEnergy <= DBL_MIN) return;
  DefineMaterial(couple);
  SetupParticle(dp->GetDefinition());
  G4double ekin = std::max(lowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  //G4cout << "G4eCoulombScatteringModel::SampleSecondaries e(MeV)= " 
  //	 << kinEnergy << "  " << particle->GetParticleName() << G4endl;
 
  // Choose nucleus
  currentElement = SelectRandomAtom(couple,particle,ekin,cutEnergy,ekin);

  SetupTarget(currentElement->GetZ(),ekin);
  
  G4double cost = SampleCosineTheta();
  G4double z1   = 1.0 - cost;
  if(z1 < 0.0) return;

  G4double sint = sqrt(z1*(1.0 + cost));
  
  //G4cout<<"## Sampled sint= " << sint << "  Z= " << targetZ 
  //	<< "  screenZ= " << screenZ << " cn= " << formfactA << G4endl;
  
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   

  // recoil sampling assuming a small recoil
  // and first order correction to primary 4-momentum
  if(lowEnergyLimit < kinEnergy) {
    G4int ia = SelectIsotopeNumber(currentElement);
    G4double Trec = z1*mom2/(amu_c2*G4double(ia));
    G4double th = 
      std::min(recoilThreshold,
	       targetZ*currentElement->GetIonisation()->GetMeanExcitationEnergy());

    if(Trec > th) {
      G4int iz = G4int(targetZ);
      G4ParticleDefinition* ion = theParticleTable->FindIon(iz, ia, 0, iz);
      Trec = z1*mom2/ion->GetPDGMass();
      if(Trec < kinEnergy) {
	G4ThreeVector dir = (direction - newDirection).unit();
	G4DynamicParticle* newdp = new G4DynamicParticle(ion, dir, Trec);
	fvect->push_back(newdp);
	fParticleChange->SetProposedKineticEnergy(kinEnergy - Trec);
      }
    }
  }
 
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::SampleCosineTheta()
{
  G4double costm = cosTetMaxNuc2;
  G4double formf = formfactA;
  G4double prob  = 0.0; 
  G4double xs = CrossSectionPerAtom();
  if(xs > 0.0) prob = elecXSection/xs;

  // scattering off e or A?
  if(G4UniformRand() < prob) {
    costm = cosTetMaxElec2;
    formf = 0.0;
  }

  /*
  G4cout << "SampleCost: e(MeV)= " << tkin 
  	 << " ctmin= " << cosThetaMin
  	 << " ctmaxN= " << cosTetMaxNuc
  	 << " ctmax= " << costm
  	 << " Z= " << targetZ << " A= " << targetA
  	 << G4endl;
  */
  if(costm >= cosTetMinNuc) return 2.0; 

  G4double x1 = 1. - cosTetMinNuc + screenZ;
  G4double x2 = 1. - costm + screenZ;
  G4double x3 = cosTetMinNuc - costm;
  G4double grej, z1; 
  do {
    z1 = x1*x2/(x1 + G4UniformRand()*x3) - screenZ;
    grej = 1.0/(1.0 + formf*z1);
  } while ( G4UniformRand() > grej*grej );  

  //G4cout << "z= " << z1 << " cross= " << nucXSection/barn 
  // << " crossE= " << elecXSection/barn << G4endl;

  return 1.0 - z1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


