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
// $Id: G4CoulombScatteringModel.cc,v 1.8 2007/05/22 17:34:36 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.10.06 V.Ivanchenko use inheritance from G4eCoulombScatteringModel
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4CoulombScatteringModel.hh"
#include "Randomize.hh"
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
  : G4eCoulombScatteringModel(thetaMin,thetaMax,build,tlim,nam)
{
  theMatManager    = G4NistManager::Instance();
  theParticleTable = G4ParticleTable::GetParticleTable();
  theProton        = G4Proton::Proton();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScatteringModel::~G4CoulombScatteringModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CoulombScatteringModel::CalculateCrossSectionPerAtom(
		             const G4ParticleDefinition* p,      
			     G4double kinEnergy, 
			     G4double Z)
{
  G4double cross= 0.0;
  G4int iz      = G4int(Z);
  G4double m    = p->GetPDGMass();
  G4double mom2 = kinEnergy*(kinEnergy + 2.0*m);
  G4double mass2= m*m;
  G4double m1   = theMatManager->GetAtomicMassAmu(iz)*amu_c2;
  G4double etot = kinEnergy + m + m1;
  G4double ptot = sqrt(mom2);
  G4double bet  = ptot/etot;
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double momCM  = gam*(ptot - bet*etot);
  G4double momCM2 = momCM*momCM;
  G4double costm = std::max(cosThetaMax, 1.0 - 0.5*q2Limit/momCM2);
  if(1 == iz && p == theProton) costm = std::max(0.0, costm);

  // Cross section in CM system 
  if(costm < cosThetaMin) {
    G4double q        = p->GetPDGCharge()/eplus;
    G4double q2       = q*q;
    G4double invbeta2 = 1.0 +  mass2/momCM2;
    G4double A = ScreeningParameter(Z, q2, momCM2, invbeta2);
    G4double a = 2.0*A + 1.0;
    G4double f = q * m1 /(m + m1);
    cross = coeff*f*f*Z*(Z + 1.0)*invbeta2*(cosThetaMin - costm)/
      ((a - cosThetaMin)*(a - costm)*momCM2);
  }
  //G4cout << "p= " << mom << " momCM= " << momCM << "  Z= " << Z << "  A= " << A 
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
      G4cout << "G4CoulombScatteringModel::SelectIsotope WARNING: "
	     << "abandance vector for"
	     << elm->GetName() << " is not normalised to unit" << G4endl;
    } else {
      N = G4double(elm->GetIsotope(idx)->GetN());
    }
  }
  return N;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScatteringModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						 const G4MaterialCutsCouple* couple,
						 const G4DynamicParticle* dp,
						 G4double,
						 G4double)
{
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();

  const G4Element* elm = 
    SelectRandomAtom(aMaterial, p, dp->GetKineticEnergy());
  G4double Z  = elm->GetZ();
  G4double N  = SelectIsotope(elm);
  G4int iz    = G4int(Z);
  G4int in    = G4int(N + 0.5);
  G4double m2 = theParticleTable->GetIonTable()->GetNucleusMass(iz, in);
  G4double m1 = dp->GetMass();

  G4double q  = p->GetPDGCharge()/eplus;
  G4double q2 = q*q;

  // Transformation to CM system
  G4LorentzVector lv1 = dp->Get4Momentum();
  G4ThreeVector dir   = dp->GetMomentumDirection();
  G4LorentzVector lv2(0.0,0.0,0.0,m2);
  G4LorentzVector lv = lv1 + lv2;
  G4ThreeVector bst  = lv.boostVector();
  lv1.boost(-bst);
  lv2.boost(-bst);
  G4ThreeVector p1   = lv1.vect();
  G4double momCM2    = p1.mag2();
  G4double invbeta2  = 1.0 + m1*m1/momCM2;
  G4double A = ScreeningParameter(Z, q2, momCM2, invbeta2);
  G4double a = 2.0*A + 1.0;

  G4double costm = std::max(cosThetaMax, 1.0 - 0.5*q2Limit/momCM2);
  if(1 == iz && p == theProton) costm = std::max(0.0, costm);
  if(costm > cosThetaMin) return; 

  G4double x   = G4UniformRand();
  G4double y   = (a + 1.0 - cosThetaMin)/(cosThetaMin - costm);
  G4double st2 = 0.5*(y*(1.0 - costm) - a*x)/(y + x); 
  if(st2 < 0.0 || st2 > 1.0) {
    G4cout << "G4CoulombScatteringModel::SampleSecondaries WARNING st2= " 
	   << st2 << G4endl;
    st2 = 0.0;
  }

  G4double tet = 2.0*asin(sqrt(st2));
  G4double cost= cos(tet);
  G4double sint= sin(tet);

  G4double phi  = twopi * G4UniformRand();

  G4ThreeVector v1(cos(phi)*sint,sin(phi)*sint,cost);
  G4double p1tot = sqrt(momCM2);
  //  v1.rotateUz(p1);
  G4LorentzVector lfv1(v1.x()*p1tot,v1.y()*p1tot,v1.z(),lv1.e());
  lfv1.boost(bst);

  G4LorentzVector lfv2 = lv - lfv1;

  G4ThreeVector newdir = lfv1.vect().unit();
  fParticleChange->ProposeMomentumDirection(newdir);   
  G4double ekin = lfv1.e() - m1;
  if(ekin < 0.0) ekin = 0.0;
  fParticleChange->SetProposedKineticEnergy(ekin);

  ekin = lfv2.e() - m2;
  if(ekin > Z*aMaterial->GetIonisation()->GetMeanExcitationEnergy()) {
    G4ParticleDefinition* ion = theParticleTable->GetIon(iz, in, 0.0);
    G4DynamicParticle* newdp  = new G4DynamicParticle(ion, lfv2);
    fvect->push_back(newdp);
  } else if(ekin > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(ekin);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


