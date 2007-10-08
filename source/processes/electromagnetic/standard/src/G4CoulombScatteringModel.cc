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
// $Id: G4CoulombScatteringModel.cc,v 1.21 2007-10-08 10:35:45 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScatteringModel::~G4CoulombScatteringModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CoulombScatteringModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy, 
				G4double Z, 
				G4double A, 
				G4double cutEnergy,
				G4double maxEnergy)
{
  // Lab system
  G4double cross= 0.0;
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  SetupTarget(Z, A, ekin);
  G4double ecross = 
    ComputeElectronXSectionPerAtom(p,tkin,Z,A,cutEnergy,maxEnergy);

  // CM system
  G4int iz      = G4int(Z);
  G4double m1   = theMatManager->GetAtomicMassAmu(iz)*amu_c2;
  G4double etot = tkin + mass;
  G4double ptot = sqrt(mom2);
  G4double bet  = ptot/(etot + m1);
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double momCM= gam*(ptot - bet*etot);

  //  G4cout << "ptot= " << ptot << " etot= " << etot << " beta= " 
  //	 << bet << " gam= " << gam << " Z= " << Z << " A= " << A << G4endl;
  // G4cout << " CM. mom= " << momCM  << " m= " << m 
  // << " m1= " << m1 << " iz= " << iz <<G4endl;

  G4double momCM2 = momCM*momCM;
  cosTetMaxNuc = std::max(cosThetaMax, 1.0 - 0.5*q2Limit/momCM2);
  if(1.5 > targetA && p == theProton && cosTetMaxNuc < 0.0) cosTetMaxNuc = 0.0;
  //G4cout << " ctmax= " << cosTetMaxNuc << " ctmin= " << cosThetaMin << G4endl;  

  // Cross section in CM system 
  if(cosTetMaxNuc < cosThetaMin) {
    G4double effmass = mass*m1/(mass + m1);
    G4double x1 = 1.0 - cosThetaMin + screenZ;
    G4double x2 = 1.0 - cosTetMaxNuc + screenZ;
    cross = coeff*Z*Z*chargeSquare*(1.0 +  effmass*effmass/momCM2)
      *(1./x1 - 1./x2 - formfactA*(2.*log(x2/x1) - 1.))/momCM2;
    //G4cout << "XS: x1= " << x1 << " x2= " << x2 << " cross= " << cross << G4endl;
    //G4cout << "momCM2= " << momCM2 << " invbeta2= " << invbeta2 
    //       << " coeff= " << coeff << G4endl;
  }
  cross += ecross;
  if(cross < 0.0) cross = 0.0;
  //  G4cout << "p= " << sqrt(mom2) << " momCM= " << momCM 
  //         << "  Z= " << Z << "  A= " << A 
  //	     << " cross= " << cross << " m1(GeV)=  " << m1/GeV <<G4endl;
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

void G4CoulombScatteringModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double cut, G4double tmax)
{
  const G4Material* aMaterial = couple->GetMaterial();
  const G4ParticleDefinition* p = dp->GetDefinition();
  G4double kinEnergy = dp->GetKineticEnergy();

  // Select isotope and setup
  SetupParticle(p);
  const G4Element* elm = SelectRandomAtom(aMaterial, p, kinEnergy);
  G4double Z  = elm->GetZ();
  G4double A  = SelectIsotope(elm);
  G4int iz    = G4int(Z);
  G4int ia    = G4int(A + 0.5);

  G4double cross = ComputeCrossSectionPerAtom(p,kinEnergy,Z,A,cut,tmax);

  G4double costm = cosTetMaxNuc;
  G4double formf = formfactA;
  if(G4UniformRand()*cross < elecXSection) {
    costm = cosTetMaxElec;
    formf = 0.0;
  }

  //  G4cout << "SampleSec: Ekin= " << kinEnergy << " m1= " << m1 
  // << " Z= "<< Z << " A= " <<A<< G4endl; 

  if(costm >= cosThetaMin) return; 

  // kinematics in CM system
  G4double m1   = theParticleTable->GetIonTable()->GetNucleusMass(iz, ia);
  G4double etot = kinEnergy + mass;
  G4double ptot = sqrt(mom2);
  G4double bet  = ptot/(etot + m1);
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double pCM  = gam*(ptot - bet*etot);
  G4double eCM  = gam*(etot - bet*ptot);

  G4double x1 = 1. - cosThetaMin + screenZ;
  G4double x2 = 1. - costm;
  G4double x3 = cosThetaMin - costm;
  G4double cost, st2, grej,  z, z1; 
  do {
    z  = G4UniformRand()*x3;
    z1 = (x1*x2 - screenZ*z)/(x1 + z);
    if(z1 < 0.0) z1 = 0.0;
    else if(z1 > 2.0) z1 = 2.0;
    cost = 1.0 - z1;
    st2  = z1*(2.0 - z1);
    grej = 1.0/(1.0 + formfactA*st2);
    //G4cout << "majorant= " << grej << " cost= " << cost << G4endl;
  } while ( G4UniformRand() > grej*grej );  

  G4double sint= sqrt(st2);

  G4double phi = twopi * G4UniformRand();

  // projectile after scattering
  G4double pzCM = pCM*cost;
  G4ThreeVector v1(pCM*cos(phi)*sint,pCM*sin(phi)*sint,gam*(pzCM + bet*eCM));
  G4ThreeVector dir = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection = v1.unit();
  newDirection.rotateUz(dir);   
  fParticleChange->ProposeMomentumDirection(newDirection);   
  G4double elab = gam*(eCM + bet*pzCM);
  G4double ekin = elab - mass;
  if(ekin < 0.0) ekin = 0.0;
  G4double plab = sqrt(ekin*(ekin + 2.0*mass));
  fParticleChange->SetProposedKineticEnergy(ekin);

  // recoil
  G4double erec = kinEnergy - ekin;
  if(erec > Z*aMaterial->GetIonisation()->GetMeanExcitationEnergy()) {
    G4ParticleDefinition* ion = theParticleTable->FindIon(iz, ia, 0, iz);
    G4ThreeVector p2 = (ptot*dir - plab*newDirection).unit();
    G4DynamicParticle* newdp  = new G4DynamicParticle(ion, p2, erec);
    fvect->push_back(newdp);
  } else if(erec > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(erec);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


