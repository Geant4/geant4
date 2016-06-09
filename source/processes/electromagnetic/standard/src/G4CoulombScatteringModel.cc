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
// $Id: G4CoulombScatteringModel.cc,v 1.2 2006/06/29 19:52:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CoulombScatteringModel
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

#include "G4CoulombScatteringModel.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4CoulombScatteringModel::G4CoulombScatteringModel(
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
  theMatManager    = G4NistManager::Instance();
  theParticleTable = G4ParticleTable::GetParticleTable();
  theProton        = G4Proton::Proton();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScatteringModel::~G4CoulombScatteringModel()
{
  if(theCrossSectionTable) {
    theCrossSectionTable->clearAndDestroy();
    delete theCrossSectionTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
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

G4double G4CoulombScatteringModel::CalculateCrossSectionPerAtom(
		             const G4ParticleDefinition* p,      
			     G4double mom2, 
			     G4double Z)
{
  G4double cross= 0.0;
  G4int iz      = G4int(Z);
  G4double m    = p->GetPDGMass();
  G4double mass2= m*m;
  G4double q    = p->GetPDGCharge()/eplus;
  G4double m1   = theMatManager->GetAtomicMassAmu(iz)*amu_c2;
  G4double etot = sqrt(mom2 + mass2) + m1;
  G4double ptot = sqrt(mom2);
  G4double bet  = ptot/etot;
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double mom  = gam*(ptot - bet*etot);
  G4double momentum2 = mom*mom;
  G4double costm = std::max(cosThetaMax, 1.0 - q2Limit/2.0*momentum2);
  if(1 == iz && p == theProton) costm = std::max(0.0, costm);

  // Cross section in CM system 
  if(costm < cosThetaMin) {
    G4double invbeta2 = 1.0 +  mass2/momentum2;
    G4double a = 2.0*pow(Z,0.666666667)*a0*
      (1.13 + 3.76*invbeta2*Z*Z*fine_structure_const*fine_structure_const)/momentum2 + 1.0;
    G4double f = q * Z * m1 /(m + m1);
    cross = coeff*f*f*(cosThetaMin - costm)/((a - cosThetaMin)*(a - costm));
  }
  //G4cout << "p= " << mom << "  Z= " << Z << "  a= " << a 
  //<< " cross= " << cross << " m1(GeV)=  " << m1/GeV <<G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CoulombScatteringModel::SelectIsotope(const G4Element* elm)
{
  G4double N = elm->GetN(); 
  G4int ni   = elm->GetNumberOfIsotopes();
  if(ni > 0) {
    G4double* ab = elm->GetRelativeAbundanceVector();
    G4double x = G4UniformRand();
    G4int idx;
    for(idx=0; idx<ni; idx++) { 
      x -= ab[idx];
      if (x <= 0.0) break;
    }
    if(idx >= ni) {
      G4cout << "G4CoulombScatteringModel::SelectIsotope WARNING: abandance vector for"
	     << elm->GetName() << " is not normalised to unit" << G4endl;
    } else {
      N = G4double(elm->GetIsotope(idx)->GetN());
    }
  }
  return N;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4CoulombScatteringModel::SampleSecondaries(
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double,
                                   G4double)
{
  std::vector<G4DynamicParticle*>* fvect = 0; 
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();

  const G4Element* elm = SelectRandomAtom(aMaterial, p, dp->GetKineticEnergy());
  G4double Z  = elm->GetZ();
  G4double N  = SelectIsotope(elm);
  G4int iz    = G4int(Z);
  G4int in    = G4int(N + 0.5);
  G4double m2 = theParticleTable->GetIonTable()->GetNucleusMass(iz, in);
  G4double m1 = dp->GetMass();

  // Transformation to CM system
  G4LorentzVector lv1 = dp->Get4Momentum();
  G4LorentzVector lv2(0.0,0.0,0.0,m2);
  G4LorentzVector lv = lv1 + lv2;
  G4ThreeVector bst  = lv.boostVector();
  lv1.boost(-bst);
  lv2.boost(-bst);
  G4ThreeVector p1   = lv1.vect();
  G4double momentum2 = p1.mag2();
  G4double invbeta2  = 1.0 + m1*m1/momentum2;

  G4double a = 2.*pow(Z,0.666666667)*a0*
    (1.13 + 3.76*invbeta2*Z*Z*fine_structure_const*fine_structure_const)/momentum2 + 1.0;

  G4double costm = std::max(cosThetaMax, 1.0 - q2Limit/2.0*momentum2);
  if(1 == iz && p == theProton) costm = std::max(0.0, costm);
  if(costm > cosThetaMin) return fvect; 

  G4double cost = a - (a - cosThetaMin)*(a - costm)/
    (a - cosThetaMin + G4UniformRand()*(cosThetaMin - costm));
  if(std::abs(cost) > 1.) {
    G4cout << "G4CoulombScatteringModel::SampleSecondaries WARNING cost= " << cost << G4endl;
    if(cost < -1.) cost = -1.0;
    else           cost =  1.0;
  }
  G4double sint = sqrt((1.0 + cost)*(1.0 - cost));

  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector v1(cos(phi)*sint,sin(phi)*sint,cost);
  G4double p1tot = sqrt(momentum2);
  v1.rotateUz(p1);
  G4LorentzVector lfv1(v1.x()*p1tot,v1.y()*p1tot,v1.z(),lv1.e());
  G4LorentzVector lfv2 = lv1 + lv2 - lfv1;

  lfv1.boost(bst);
  lfv2.boost(bst);

  fParticleChange->ProposeMomentumDirection(lfv1.vect().unit());   
  G4double ekin = lfv1.e() - m1;
  if(ekin < 0.0) ekin = 0.0;
  fParticleChange->SetProposedKineticEnergy(ekin);

  ekin = lfv2.e() - m2;
  if(ekin > Z*aMaterial->GetIonisation()->GetMeanExcitationEnergy()) {
    fvect = new std::vector<G4DynamicParticle*>;
    G4ParticleDefinition* ion = theParticleTable->GetIon(iz, in, 0.0);
    G4DynamicParticle* newdp  = new G4DynamicParticle(ion, lfv2);
    fvect->push_back(newdp);
  } 
 
  return fvect;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


