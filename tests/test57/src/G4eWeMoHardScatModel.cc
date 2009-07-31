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
// $Id: G4eWeMoHardScatModel.cc,v 1.1 2009-07-31 15:31:32 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eWeMoHardScatModel
//
// Author:        V. Grichine based on G4eCoulombScatteringModel 
//
// Creation date: 31.07.2009
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

#include "G4eWeMoHardScatModel.hh"
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

G4double G4eWeMoHardScatModel::fScreenRSquare[] = {0.0};
G4double G4eWeMoHardScatModel::fFormFactor[]    = {0.0};

using namespace std;

G4eWeMoHardScatModel::G4eWeMoHardScatModel(const G4String& nam)
  : G4VEmModel(nam),
    fCosThetaMin(1.0),
    fCosThetaMax(-1.0),
    fq2Limit(TeV*TeV),
    fAlpha2(fine_structure_const*fine_structure_const),
    fFacLim(100.0),
    fInitialised(false)
{
  fNistManager = G4NistManager::Instance();
  theParticleTable = G4ParticleTable::GetParticleTable();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  theProton   = G4Proton::Proton();
  fCurrentMaterial = 0; 
  fCurrentElement  = 0;
  fLowEnergyLimit = keV;
  G4double p0 = electron_mass_c2*classic_electr_radius;
  fCoeff  = twopi*p0*p0;
  fTkin = fTargetZ = fMom2 = DBL_MIN;
  fElecXSection = fNucXSection = 0.0;
  fRecoilThreshold = 100.*keV;
  feCut = DBL_MAX;
  fParticle = 0;
  fCurrentCouple = 0;

  // Thomas-Fermi screening radii
  // Formfactors from A.V. Butkevich et al., NIM A 488 (2002) 282

  if(0.0 == fScreenRSquare[0]) {
    G4double a0 = electron_mass_c2/0.88534; 
    G4double constn = 6.937e-6/(MeV*MeV);

    fScreenRSquare[0] = fAlpha2*a0*a0;
    for(G4int j=1; j<100; j++) {
      G4double x = a0*fNistManager->GetZ13(j);
      fScreenRSquare[j] = fAlpha2*x*x;
      x = fNistManager->GetA27(j); 
      fFormFactor[j] = constn*x*x;
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eWeMoHardScatModel::~G4eWeMoHardScatModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eWeMoHardScatModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector& cuts)
{
  SetupParticle(p);
  fCurrentCouple = 0;
  fElecXSection = fNucXSection = 0.0;
  fTkin = fTargetZ = fMom2 = DBL_MIN;
  feCut = feTag = DBL_MAX;
  fCosThetaMin = cos(PolarAngleLimit());
  fCurrentCuts = &cuts;
  //G4cout << "!!! G4eWeMoHardScatModel::Initialise for " 
  //	 << p->GetParticleName() << "  cos(TetMin)= " << fCosThetaMin 
  //	 << "  cos(TetMax)= " << fCosThetaMax <<G4endl;
  if(!fInitialised) {
    fInitialised = true;
    fParticleChange = GetParticleChangeForGamma();
 }
  if(fMass < GeV && fParticle->GetParticleType() != "nucleus") {
    InitialiseElementSelectors(p,cuts);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4eWeMoHardScatModel::ComputeMaxElectronScattering(G4double cutEnergy)
{
  feCut = cutEnergy;
  G4double tmax = fTkin;
  fCosTetMaxElec = 1.0;
  if(fMass > MeV) {
    G4double ratio = electron_mass_c2/fMass;
    G4double tau = fTkin/fMass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
    fCosTetMaxElec = 1.0 - std::min(cutEnergy, tmax)*electron_mass_c2/fMom2;
  } else {

    if(fParticle == theElectron) tmax *= 0.5;
    G4double t = std::min(cutEnergy, tmax);
    G4double mom21 = t*(t + 2.0*electron_mass_c2);
    G4double t1 = fTkin - t;
    //G4cout << "fTkin= " << fTkin << " t= " << t << " t1= " << t1 << G4endl;
    if(t1 > 0.0) {
      G4double mom22 = t1*(t1 + 2.0*fMass);
      G4double ctm = (fMom2 + mom22 - mom21)*0.5/sqrt(fMom2*mom22);
      //G4cout << "ctm= " << ctm << G4endl;
      if(ctm <  1.0) fCosTetMaxElec = ctm;
      if(ctm < -1.0) fCosTetMaxElec = -1.0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eWeMoHardScatModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double,
		G4double cutEnergy, G4double)
{
  //G4cout << "### G4eWeMoHardScatModel::ComputeCrossSectionPerAtom  for " 
  //  << p->GetParticleName()<<" Z= "<<Z<<" e(MeV)= "<< kinEnergy/MeV << G4endl; 
  G4double xsec = 0.0;
  SetupParticle(p);
  G4double ekin = std::max(fLowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  if(fCosTetMaxNuc < fCosTetMinNuc) {
    SetupTarget(Z, ekin);
    xsec = CrossSectionPerAtom();  
  }
  /*
  G4cout << "e(MeV)= " << ekin/MeV << "fCosTetMinNuc= " << fCosTetMinNuc
	 << " fCosTetMaxNuc= " << fCosTetMaxNuc
	 << " fCosTetMaxElec= " << fCosTetMaxElec
	 << " fScreenZ= " << fScreenZ
	 << " fFormFactA= " << fFormFactA
	 << " fCosTetMaxHad= " << fCosTetMaxHad << G4endl;
  */
  return xsec;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eWeMoHardScatModel::CrossSectionPerAtom()
{
  // This method needs initialisation before be called
  //  G4double fac = fCoeff*fTargetZ*fChargeSquare*kinFactor;

  G4double gamma2=1./(1.-1./fInvBeta2);
  G4double mu_c2 =(fMass*fTargetMass)/(fMass+fTargetMass);
  G4double fac = fCoeff*fTargetZ*fChargeSquare*fInvBeta2*fInvBeta2/(gamma2*mu_c2*mu_c2);

  fElecXSection = 0.0;
  fNucXSection  = 0.0;

  G4double x  = 1.0 - fCosTetMinNuc;
  G4double x1 = x + fScreenZ;

  if(fCosTetMaxElec2 < fCosTetMinNuc) {
    fElecXSection = fac*(fCosTetMinNuc - fCosTetMaxElec2)/
      (x1*(1.0 - fCosTetMaxElec2 + fScreenZ));
    fNucXSection  = fElecXSection;
  }

  //G4cout << "XS fTkin(MeV)= " << fTkin<<" xs= " <<fNucXSection 
  //	 << " fCostmax= " << fCosTetMaxNuc2 
  //	 << " fCostmin= " << fCosTetMinNuc << "  Z= " << fTargetZ <<G4endl;
  if(fCosTetMaxNuc2 < fCosTetMinNuc) {
    G4double s  = fScreenZ*fFormFactA;
    G4double z1 = 1.0 - fCosTetMaxNuc2 + fScreenZ;
    G4double s1 = 1.0 - s;
    G4double d  = s1/fFormFactA;
    //G4cout <<"x1= "<<x1<<" z1= " <<z1<<" s= "<<s << " d= " <<d <<G4endl;
    if(d < 0.2*x1) {
      G4double x2 = x1*x1;
      G4double z2 = z1*z1;
      x = (1.0/(x1*x2) - 1.0/(z1*z2) - d*1.5*(1.0/(x2*x2) - 1.0/(z2*z2)))/
	(3.0*fFormFactA*fFormFactA);
    } else {
      G4double x2 = x1 + d;
      G4double z2 = z1 + d;
      x = (1.0/x1 - 1.0/z1 + 1.0/x2 - 1.0/z2 - 2.0*log(z1*x2/(z2*x1))/d)/(s1*s1);
    }
    fNucXSection += fac*fTargetZ*x;
  }
  //G4cout<<" cross(bn)= "<<fNucXSection/barn<<" xsElec(bn)= "<<fElecXSection/barn
  //	<< " Asc= " << fScreenZ << G4endl; 
  
  return fNucXSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eWeMoHardScatModel::SampleSecondaries(
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
  G4double ekin = std::max(fLowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  //G4cout << "G4eWeMoHardScatModel::SampleSecondaries e(MeV)= " 
  //	 << kinEnergy << "  " << fParticle->GetParticleName() << G4endl;
 
  // Choose nucleus
  fCurrentElement = SelectRandomAtom(couple,fParticle,ekin,cutEnergy,ekin);

  SetupTarget(fCurrentElement->GetZ(),ekin);
  
  G4double cost = SampleCosineTheta();
  G4double z1   = 1.0 - cost;
  if(z1 < 0.0) return;

  G4double sint = sqrt(z1*(1.0 + cost));
  
  //G4cout<<"## Sampled sint= " << sint << "  Z= " << fTargetZ 
  //	<< "  fScreenZ= " << fScreenZ << " cn= " << fFormFactA << G4endl;
  
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   

  // recoil sampling assuming a small recoil
  // and first order correction to primary 4-momentum
  if(fLowEnergyLimit < kinEnergy) {
    G4int ia = SelectIsotopeNumber(fCurrentElement);
    G4double Trec = z1*fMom2/(amu_c2*G4double(ia));

    if(Trec > fRecoilThreshold) 
    {
      G4ParticleDefinition* ion = theParticleTable->FindIon(fiz, ia, 0, fiz);
      Trec = z1*fMom2/ion->GetPDGMass();
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

G4double G4eWeMoHardScatModel::SampleCosineTheta()
{
  G4double costm = fCosTetMaxNuc2;
  G4double formf = fFormFactA;
  G4double prob  = 0.0; 
  G4double xs = CrossSectionPerAtom();
  if(xs > 0.0) prob = fElecXSection/xs;

  // scattering off e or A?
  if(G4UniformRand() < prob) {
    costm = fCosTetMaxElec2;
    formf = 0.0;
  }

  /*
  G4cout << "SampleCost: e(MeV)= " << fTkin 
  	 << " ctmin= " << fCosThetaMin
  	 << " ctmaxN= " << fCosTetMaxNuc
  	 << " ctmax= " << costm
  	 << " Z= " << fTargetZ << " A= " << targetA
  	 << G4endl;
  */
  if(costm >= fCosTetMinNuc) return 2.0; 

  G4double x1 = 1. - fCosTetMinNuc + fScreenZ;
  G4double x2 = 1. - costm + fScreenZ;
  G4double x3 = fCosTetMinNuc - costm;
  G4double grej, z1; 
  do {
    z1 = x1*x2/(x1 + G4UniformRand()*x3) - fScreenZ;
    grej = 1.0/(1.0 + formf*z1);
  } while ( G4UniformRand() > grej*grej );  

  //G4cout << "z= " << z1 << " cross= " << fNucXSection/barn 
  // << " crossE= " << fElecXSection/barn << G4endl;

  return 1.0 - z1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


