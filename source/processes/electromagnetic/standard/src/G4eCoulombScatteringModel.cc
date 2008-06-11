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
// $Id: G4eCoulombScatteringModel.cc,v 1.52 2008-06-11 08:51:00 vnivanch Exp $
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

G4eCoulombScatteringModel::G4eCoulombScatteringModel(
  G4double thetaMin, G4double thetaMax, G4double tlim, const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(1.0),
    cosThetaMax(-1.0),
    q2Limit(tlim),
    lowKEnergy(keV),
    highKEnergy(TeV),
    alpha2(fine_structure_const*fine_structure_const),
    faclim(100.0),
    isInitialised(false)
{
  if(thetaMin > 0.0) cosThetaMin = cos(thetaMin);
  if(thetaMax < pi)  cosThetaMax = cos(thetaMax);
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
  tkin = targetZ = targetA = mom2 = DBL_MIN;
  elecXSection = nucXSection = 0.0;
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
  tkin = targetZ = targetA = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  //  G4cout << "!!! G4eCoulombScatteringModel::Initialise for " 
  // << p->GetParticleName() << "  cos(TetMin)= " << cosThetaMin 
  // << "  cos(TetMax)= " << cosThetaMax <<G4endl;
  if(!isInitialised) {
    isInitialised = true;

    if(pParticleChange)
      fParticleChange = 
	reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChange = new G4ParticleChangeForGamma();
  } else {
    return;
  }

  currentCuts = &cuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4eCoulombScatteringModel::ComputeMaxElectronScattering(G4double cutEnergy)
{
  ecut = cutEnergy;
  G4double tmax = tkin;
  if(particle == theElectron) tmax *= 0.5;
  else if(particle != thePositron) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = tkin/mass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
  }
  cosTetMaxElec = 1.0;
  G4double t = std::min(cutEnergy, tmax);
  G4double mom21 = t*(t + 2.0*electron_mass_c2);
  G4double t1 = tkin - t;
  if(t1 > 0.0) {
    G4double mom22 = t1*(t1 + 2.0*mass);
    G4double ctm = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
    if(ctm < 1.0) cosTetMaxElec = ctm;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double A,
		G4double cutEnergy, G4double)
{
  //  if(p == particle && kinEnergy == tkin && Z == targetZ &&
  //   A == targetA && cutEnergy == ecut) return nucXSection;

  //G4cout << "### G4eCoulombScatteringModel::ComputeCrossSectionPerAtom  for " 
  // << p->GetParticleName() << " Z= " << Z << " A= " << A 
  // << " e= " << kinEnergy << G4endl; 
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  SetupTarget(Z, A, ekin);

  return CrossSectionPerAtom();
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

  if(cosTetMaxElec < cosTetMinNuc) {
    elecXSection = fac*(cosTetMinNuc - cosTetMaxElec)/
      (x1*(1.0 - cosTetMaxElec + screenZ));
    nucXSection  = elecXSection;
  }

  G4double costm = std::max(cosTetLimit,cosTetMaxHad);
  if(costm < cosTetMinNuc) {
    G4double s  = screenZ*formfactA;
    G4double z1 = 1.0 - costm + screenZ;
    G4double d  = (1.0 - s)/formfactA;
    G4double x2 = x1 + d;
    G4double z2 = z1 + d;
    nucXSection += fac*targetZ*(1.0 - 2.0*s)*
      ((cosTetMinNuc - costm)*(1.0/(x1*z1) + 1.0/(x2*z2)) - 
       2.0*log(z1*x2/(z2*x1))/d);
  }

  //G4cout<<" cross(bn)= "<<nucXSection/barn<<" Asc= "<<screenZ<<G4endl; 
  
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
  G4double ekin = std::max(keV, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  //G4cout << "G4eCoulombScatteringModel::SampleSecondaries e(MeV)= " 
  // << kinEnergy <<G4endl;
  SelectAtomRandomly();

  SetupTarget(currentElement->GetZ(),currentElement->GetN(),tkin);
  
  G4double cost = SampleCosineTheta();
  G4double z1   = 1.0 - cost;
  if(z1 < 0.0) return;

  G4double sint = sqrt(z1*(1.0 + cost));
  
  //G4cout<<"## Sampled sint= " << sint << "  Z= " << targetZ 
  //<< "  screenZ= " << screenZ << " cn= " << formfactA << G4endl;
  
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   

  // recoil sampling assuming a small recoil
  // and first order correction to primary 4-momentum

  G4double Trec = z1*mom2/(amu_c2*targetA);
  if(Trec > targetZ*currentMaterial->GetIonisation()->GetMeanExcitationEnergy())
    {
      G4double A  = SelectIsotope();
      G4int iz    = G4int(targetZ);
      G4int ia    = G4int(A + 0.5);
      G4ParticleDefinition* ion = theParticleTable->FindIon(iz, ia, 0, iz);
      G4ThreeVector dir = (direction - newDirection).unit();
      Trec = z1*mom2/ion->GetPDGMass();
      G4DynamicParticle* newdp = new G4DynamicParticle(ion, dir, Trec);
      fvect->push_back(newdp);

      fParticleChange->SetProposedKineticEnergy(kinEnergy - Trec);
    }
 
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::SampleCosineTheta()
{
  G4double costm = cosTetLimit;
  G4double formf = formfactA;

  // scattering off e or A?
  if(G4UniformRand()*xsecn[idxelm] < xsece[idxelm]) {
    costm = cosTetMaxElec;
    formf = 0.0;
  } else if (cosTetLimit < cosTetMaxHad) {
    costm = cosTetMaxHad;
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

void G4eCoulombScatteringModel::SelectAtomRandomly()
{
  const G4ElementVector* theElementVector = 
    currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  G4int nelm = currentMaterial->GetNumberOfElements();
  G4double cut = (*currentCuts)[currentMaterialIndex];

  G4double cross = 0.0;

  G4int i = 0;
  for (; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double den = theAtomNumDensityVector[i];
    ComputeCrossSectionPerAtom(particle,tkin,elm->GetZ(),elm->GetN(),cut,tkin);
    xsece[i] = elecXSection*den;
    xsecn[i] = nucXSection*den;
    // G4cout << "i= " << i << " den= " << den << " eXS= " << elecXSection 
    //   << " nXS= " << nucXSection << G4endl;
    cross   += xsecn[i];
    xsect[i] = cross;
  }

  idxelm = 0;
  if(nelm > 1) {
    G4double qsec = cross*G4UniformRand();
    for (idxelm=0; idxelm<nelm; idxelm++) {
      if(qsec <= xsect[idxelm]) break;
    }
    if(idxelm >= nelm) idxelm = nelm - 1;
  } 

  currentElement = (*theElementVector)[idxelm];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eCoulombScatteringModel::SelectIsotope()
{
  G4double N = currentElement->GetN(); 
  G4int ni   = currentElement->GetNumberOfIsotopes();
  if(ni > 0) {
    G4double* ab = currentElement->GetRelativeAbundanceVector();
    G4double x = G4UniformRand();
    G4int idx;
    for(idx=0; idx<ni; idx++) { 
      x -= ab[idx];
      if (x <= 0.0) break;
    }
    if(idx >= ni) {
      G4cout << "G4CoulombScatteringModel::SelectIsotope WARNING: "
	     << "abandance vector for"
	     << currentElement->GetName() << " is not normalised to unit" << G4endl;
    } else {
      N = G4double(currentElement->GetIsotope(idx)->GetN());
    }
  }
  return N;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


