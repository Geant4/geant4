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
// $Id: G4eCoulombScatteringModel.cc,v 1.18 2007-08-15 09:24:55 vnivanch Exp $
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
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  theProton   = G4Proton::Proton();
  a0 = alpha2*electron_mass_c2*electron_mass_c2/(0.885*0.885);
  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff  = twopi*p0*p0;
  constn = 6.937e-6/(MeV*MeV);
  tkin = targetZ = targetA = mom2 = DBL_MIN;
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
  //  G4cout << "!!! G4eCoulombScatteringModel::Initialise" << G4endl;
  if(!isInitialised) {
    isInitialised = true;

    if(pParticleChange)
      fParticleChange = 
	reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChange = new G4ParticleChangeForGamma();
  }

  // Initialise mass and charge
  if(p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    G4double q = particle->GetPDGCharge()/eplus;
    chargeSquare = q*q;
  }

  if(!buildTable || p->GetParticleName() == "GenericIon") return;

  // Compute log cross section table
  if(!theCrossSectionTable)  theCrossSectionTable = new G4PhysicsTable();

  G4PhysicsLogVector* ptrVector;
  G4double e, value;
  nbins = 2*G4int(log10(highKEnergy/lowKEnergy));

  const  G4ElementTable* elmt = G4Element::GetElementTable();
  size_t nelm =  G4Element::GetNumberOfElements();

  //  G4cout << "### G4eCoulombScatteringModel: Build table for " 
  //	 << nelm << " elements for " << p->GetParticleName() << G4endl;

  for(size_t j=0; j<nelm; j++) { 

    const G4Element* elm = (*elmt)[j]; 
    G4double Z =  elm->GetZ();
    G4double A =  elm->GetN();
    G4int iz   = G4int(Z);

    // fill table for the given Z only once
    if(index[iz] == -1) {
      ptrVector  = new G4PhysicsLogVector(lowKEnergy, highKEnergy, nbins);
      index[iz] = j;
      for(G4int i=0; i<=nbins; i++) {
	e     = ptrVector->GetLowEdgeEnergy( i ) ;
	value = CalculateCrossSectionPerAtom(p, e, Z, A);  
	ptrVector->PutValue( i, log(value) );
      }
      theCrossSectionTable->insert(ptrVector);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eCoulombScatteringModel::CalculateCrossSectionPerAtom(
		             const G4ParticleDefinition* p,      
			     G4double kinEnergy, 
			     G4double Z, G4double A)
{
  G4double cross = 0.0;
  SetupTarget(Z, A, std::max(keV, kinEnergy));

  if(cosTetMaxNuc < cosThetaMin) {
    G4double x1 = 1.0 - cosThetaMin  + screenZ;
    G4double x2 = 1.0 - cosTetMaxNuc + screenZ;
    cross = coeff*Z*(Z + 1.)*chargeSquare*invbeta2
      *(1./x1 - 1./x2 - formfactA*(2.*log(x2/x1) - 1.))/mom2;
  }
  //  G4cout << "CalculateCrossSectionPerAtom: e(MeV)= " << tkin 
  //	 << " cross(b)= " << cross/barn
  //	 << " ctmin= " << cosThetaMin
  //	 << " ctmax= " << cosTetMaxNuc       
  //	 << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eCoulombScatteringModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
						  const G4MaterialCutsCouple* couple,
						  const G4DynamicParticle* dp,
						  G4double,
						  G4double)
{
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();
  
  G4double kinEnergy = dp->GetKineticEnergy();

  const G4Element* elm = SelectRandomAtom(aMaterial, p, kinEnergy);
  G4double Z  = elm->GetZ();
  G4double A  = elm->GetN();
  SetupTarget(Z, A, kinEnergy);

  //  G4cout << "G4eCoulombScatteringModel::SampleSecondaries: e(MeV)= " << tkin 
  //	 << " ctmin= " << cosThetaMin
  //	 << " ctmax= " << cosTetMaxNuc 
  //	 << " Z= " << Z << " A= " << A
  //	 << G4endl;

  if(cosTetMaxNuc >= cosThetaMin) return; 

  G4double x1 = 1. - cosThetaMin + screenZ;
  G4double x2 = 1. - cosTetMaxNuc;
  G4double x3 = cosThetaMin - cosTetMaxNuc;
  G4double cost, st2, grej,  z, z1; 
  do {
    z  = G4UniformRand()*x3;
    z1 = (x1*x2 - screenZ*z)/(x1 + z);
    if(z1 < 0.0) z1 = 0.0;
    else if(z1 > 2.0) z1 = 2.0;
    cost = 1.0 - z1;
    st2  = z1*(2.0 - z1);
    grej = 1.0/(1.0 + formfactA*st2);
  } while ( G4UniformRand() > grej*grej );  
  
  G4double sint= sqrt(st2);

  //  G4cout<<"SampleSecondaries: e(MeV)= " << kinEnergy
  //	<< " cost= " << cost << "  Z= " << Z << "  a= " << screenZ 
  //	<< " cn= " << formfactA
  //	<< G4endl;

  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(newDirection);   
 
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


