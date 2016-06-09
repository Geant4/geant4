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
// $Id: G4eCoulombScatteringModel.cc,v 1.40 2008/01/07 08:32:01 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01-patch-01 $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eCoulombScatteringModel::G4eCoulombScatteringModel(
  G4double thetaMin, G4double thetaMax, G4bool build, 
  G4double tlim, const G4String& nam)
  : G4VEmModel(nam),
    cosThetaMin(cos(thetaMin)),
    cosThetaMax(cos(thetaMax)),
    q2Limit(tlim),
    theCrossSectionTable(0),
    lowKEnergy(keV),
    highKEnergy(TeV),
    alpha2(fine_structure_const*fine_structure_const),
    faclim(100.0),
    nbins(12),
    nmax(100),
    buildTable(build),
    isInitialised(false)
{
  fNistManager = G4NistManager::Instance();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  theProton   = G4Proton::Proton();
  a0 = alpha2*electron_mass_c2*electron_mass_c2/(0.885*0.885);
  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff  = twopi*p0*p0;
  constn = 6.937e-6/(MeV*MeV);
  tkin = targetZ = targetA = mom2 = DBL_MIN;
  elecXSection = nucXSection = 0.0;
  ecut = DBL_MAX;
  particle = 0;
  for(size_t j=0; j<100; j++) {index[j] = -1;} 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eCoulombScatteringModel::~G4eCoulombScatteringModel()
{
  if(theCrossSectionTable) {
    theCrossSectionTable->clearAndDestroy();
    delete theCrossSectionTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
					   const G4DataVector&)
{
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

  if(p->GetParticleType() == "nucleus") buildTable = false;
  if(!buildTable) return;

  // Compute log cross section table per atom
  if(!theCrossSectionTable) theCrossSectionTable = new G4PhysicsTable();
  
  nbins = 2*G4int(log10(highKEnergy/lowKEnergy));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double A,
		G4double cutEnergy, G4double)
{
  if(p == particle && kinEnergy == tkin && Z == targetZ &&
     A == targetA && cutEnergy == ecut) return nucXSection;

  //G4cout << "### G4eCoulombScatteringModel::ComputeCrossSectionPerAtom  for " 
  //	 << p->GetParticleName() << " Z= " << Z << " A= " << A 
  //	 << " e= " << kinEnergy << G4endl; 

  nucXSection = ComputeElectronXSectionPerAtom(p,kinEnergy,Z,A,cutEnergy);

  // nuclear cross section
  if(theCrossSectionTable) {
    G4bool b;
    G4int iz  = G4int(Z);
    G4int idx = index[iz];

    // compute table for given Z
    if(-1 == idx) {
      idx = theCrossSectionTable->size();
      index[iz] = idx;
      G4PhysicsLogVector* ptrVector
	= new G4PhysicsLogVector(lowKEnergy, highKEnergy, nbins);
      //  G4cout << "New vector Z= " << iz << " A= " << A << " idx= " << idx << G4endl; 
      G4double e, value;
      for(G4int i=0; i<=nbins; i++) {
	e     = ptrVector->GetLowEdgeEnergy( i ) ;
	value = CalculateCrossSectionPerAtom(p, e, Z, A);  
	ptrVector->PutValue( i, log(value) );
      }
      theCrossSectionTable->push_back(ptrVector);
    }

      // take value from the table
    nucXSection += 
      std::exp((((*theCrossSectionTable)[idx]))->GetValue(kinEnergy, b));

    // compute value from scratch
  } else nucXSection += CalculateCrossSectionPerAtom(p, kinEnergy, Z, A);
  
  //  G4cout << " cross(bn)= " << nucXSection/barn << G4endl; 
  
  if(nucXSection < 0.0) nucXSection = 0.0;
  return nucXSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eCoulombScatteringModel::ComputeElectronXSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z,
		G4double A,
		G4double cutEnergy)
{
  if(p == particle && kinEnergy == tkin && Z == targetZ &&
     cutEnergy == ecut) return elecXSection;
  ecut = cutEnergy;
  elecXSection = 0.0;
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  //G4double ekin = kinEnergy;
  SetupTarget(Z, A, ekin);

  G4double tmax = tkin;
  if(p == theElectron) tmax *= 0.5;
  else if(p != thePositron) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = tkin/mass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
  }

  cosTetMaxElec = cosTetMaxNuc;
  G4double t = std::min(cutEnergy, tmax);
  G4double mom21 = t*(t + 2.0*electron_mass_c2);
  G4double t1 = tkin - t;
  if(t1 > 0.0) {
    G4double mom22 = t1*(t1 + 2.0*mass);
    G4double ctm = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
    if(ctm > cosTetMaxElec && ctm <= 1.0) cosTetMaxElec = ctm;
  }

  if(cosTetMaxElec < cosThetaMin) {
    G4double x1 = 1.0 - cosThetaMin  + screenZ;
    G4double x2 = 1.0 - cosTetMaxElec + screenZ;
    elecXSection = coeff*Z*chargeSquare*invbeta2*
      (cosThetaMin - cosTetMaxElec)/(x1*x2*mom2);
  }
  //  G4cout << "cut= " << ecut << " e= " << tkin 
  // << " croosE(barn)= " << elecXSection/barn 
  // << " cosEl= " << cosTetMaxElec << " costmin= " << cosThetaMin << G4endl;
  return elecXSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eCoulombScatteringModel::CalculateCrossSectionPerAtom(
		             const G4ParticleDefinition* p,
			     G4double kinEnergy, 
			     G4double Z, G4double A)
{
  G4double cross = 0.0;
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  //G4double ekin = kinEnergy;
  SetupTarget(Z, A, ekin);

  if(cosTetMaxNuc < cosThetaMin) {
    G4double x1 = 1.0 - cosThetaMin;
    G4double x2 = 1.0 - cosTetMaxNuc;
    G4double x3 = cosThetaMin - cosTetMaxNuc;
    G4double z1 = x1 + screenZ;
    G4double z2 = x2 + screenZ;
    G4double d  = 1.0/formfactA - screenZ;
    G4double d1 = 1.0 - formfactA*screenZ;
    G4double zn1= x1 + d;
    G4double zn2= x2 + d;
    cross = coeff*Z*Z*chargeSquare*invbeta2
      *(x3/(z1*z2) + x3/(zn1*zn2) + 
	2.0*std::log(z1*zn2/(z2*zn1))/d) / (mom2*d1*d1);
  }
  
  //  G4cout << "CalculateCrossSectionPerAtom: e(MeV)= " << tkin 
  // << " cross(b)= " << cross/barn << " ctmin= " << cosThetaMin
  // << " ctmax= " << cosTetMaxNuc << G4endl;
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::SampleSecondaries(
                std::vector<G4DynamicParticle*>*,
		const G4MaterialCutsCouple* couple,
		const G4DynamicParticle* dp,
		G4double cutEnergy,
		G4double maxEnergy)
{
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();
  G4double kinEnergy = dp->GetKineticEnergy();

  // Select atom and setup
  SetupParticle(p);
  const G4Element* elm = 
    SelectRandomAtom(aMaterial,p,kinEnergy,cutEnergy,maxEnergy);
  G4double Z  = elm->GetZ();
  G4double A  = elm->GetN();

  G4double cross = 
    ComputeCrossSectionPerAtom(p,kinEnergy,Z,A,cutEnergy,maxEnergy);

  G4double costm = cosTetMaxNuc;
  G4double formf = formfactA;
  if(G4UniformRand()*cross < elecXSection) {
    costm = cosTetMaxElec;
    formf = 0.0;
  }
  /*
  G4cout << "G4eCoul...SampleSecondaries: e(MeV)= " << tkin 
  	 << " ctmin= " << cosThetaMin
  	 << " ctmaxN= " << cosTetMaxNuc
  	 << " ctmax= " << costm
  	 << " Z= " << Z << " A= " << A
  	 << " cross= " << cross/barn << " crossE= " << elecXSection/barn
  	 << G4endl;
  */
  if(costm >= cosThetaMin) return; 

  G4double x1 = 1. - cosThetaMin + screenZ;
  G4double x2 = 1. - costm;
  G4double x3 = cosThetaMin - costm;
  G4double grej,  z, z1; 
  do {
    z  = G4UniformRand()*x3;
    z1 = (x1*x2 - screenZ*z)/(x1 + z);
    if(z1 < 0.0) z1 = 0.0;
    else if(z1 > 2.0) z1 = 2.0;
    grej = 1.0/(1.0 + formf*z1);
  } while ( G4UniformRand() > grej*grej );  
  
  G4double cost = 1.0 - z1;
  G4double sint= sqrt(z1*(2.0 - z1));
  /*
  if(sint > 0.1) 
    G4cout<<"## SampleSecondaries: e(MeV)= " << kinEnergy
	  << " sint= " << sint << "  Z= " << Z << "  screenZ= " << screenZ 
	  << " cn= " << formf
	  << G4endl;
  */
  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   
 
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


