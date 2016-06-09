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
// $Id: G4eCoulombScatteringModel.cc,v 1.2 2006/06/29 19:53:49 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eCoulombScatteringModel.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eCoulombScatteringModel::G4eCoulombScatteringModel(
  G4double thetaMin, G4double thetaMax, G4bool build, 
  G4double tlim, const G4String& nam)
  : G4VEmModel(nam),
    theCrossSectionTable(0),
    cosThetaMin(cos(thetaMin)),
    cosThetaMax(cos(thetaMax)),
    lowMomentum(keV),
    highMomentum(MeV),
    q2Limit(tlim),
    nbins(12),
    nmax(100),
    buildTable(build),
    isInitialised(false)
{
  G4double p0 = hbarc/(Bohr_radius*0.885);
  a0 = 0.25*p0*p0;
  p0 = electron_mass_c2*classic_electr_radius;
  coeff = twopi*p0*p0;
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
  if(isInitialised) return;
  isInitialised = true;

  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForGamma();

  if(!buildTable || p->GetParticleName() == "GenericIon") return;

  // Compute cross section multiplied by Ptot^2*beta^2
  G4double mass  = p->GetPDGMass();
  G4double mass2 = mass*mass;

  theCrossSectionTable = new G4PhysicsTable(nmax);
  G4PhysicsLogVector* ptrVector;
  G4double mom2, value;
  G4double pmin = lowMomentum*lowMomentum;
  G4double pmax = highMomentum*highMomentum;
  nbins = G4int(log10(pmax/pmin)/2.0) + 1;

  for(G4int j=1; j<nmax; j++) { 

    ptrVector  = new G4PhysicsLogVector(pmin, pmax, nbins);
 
    for(G4int i=0; i<=nbins; i++) {
      mom2   = ptrVector->GetLowEdgeEnergy( i ) ;
      value  = CalculateCrossSectionPerAtom(p, mom2, j);  
      value *= mom2*mom2/(mom2 + mass2);
      ptrVector->PutValue( i, value );
    }

    theCrossSectionTable->insert(ptrVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eCoulombScatteringModel::CalculateCrossSectionPerAtom(
		             const G4ParticleDefinition* p,      
			     G4double momentum2, 
			     G4double Z)
{
  G4double cross = 0.0;
  G4double m     = p->GetPDGMass();
  G4double q     = p->GetPDGCharge()/eplus;
  G4double mass2 = m*m;
  G4double costm = std::max(cosThetaMax, 1.0 - q2Limit/2.0*momentum2);
  if(costm < cosThetaMin) {
    G4double invbeta2 = 1.0 +  mass2/momentum2;
    G4double a = 2.0*pow(Z,0.666666667)*a0*
      (1.13 + 3.76*invbeta2*Z*Z*fine_structure_const*fine_structure_const)/momentum2 + 1.0;
    cross = coeff*q*q*Z*Z*(cosThetaMin - costm)/((a - cosThetaMin)*(a - costm));
  }
  //G4cout << "p= " << sqrt(momentum2) << "  Z= " << Z << "  a= " << a 
  //<< " cross= " << cross << " " <<G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4eCoulombScatteringModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double,
                                   G4double)
{
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();
  
  G4double mass      = dp->GetMass();
  G4double kinEnergy = dp->GetKineticEnergy();
  G4double momentum2 = kinEnergy*(kinEnergy + 2.0*mass);
  G4double invbeta2  = (kinEnergy + mass)*(kinEnergy + mass)/momentum2;

  const G4Element* elm = SelectRandomAtom(aMaterial, p, kinEnergy);
  G4double Z  = elm->GetZ();

  G4double a = 2.*pow(Z,0.666666667)*a0*
    (1.13 + 3.76*invbeta2*Z*Z*fine_structure_const*fine_structure_const)/momentum2 + 1.0;
  G4double costm = std::max(cosThetaMax, 1.0 - q2Limit/2.0*momentum2);
  if(costm > cosThetaMin) return 0; 

  G4double cost = a - (a - cosThetaMin)*(a - costm)/
    (a - cosThetaMin + G4UniformRand()*(cosThetaMin - costm));
  if(std::abs(cost) > 1.) {
    G4cout << "G4eCoulombScatteringModel::SampleSecondaries WARNING cost= " << cost << G4endl;
    if(cost < -1.) cost = -1.0;
    else           cost =  1.0;
  }
  G4double sint = sqrt((1.0 + cost)*(1.0 - cost));

  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector direction = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection(cos(phi)*sint,sin(phi)*sint,cost);
  newDirection.rotateUz(direction);   

  fParticleChange->ProposeMomentumDirection(direction);   
 
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


