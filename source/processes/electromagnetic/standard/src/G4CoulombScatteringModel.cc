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
// $Id: G4CoulombScatteringModel.cc,v 1.44.2.1 2010/04/06 09:14:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-03-patch-02 $
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
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.10.06 V.Ivanchenko use inheritance from G4eCoulombScatteringModel
// 09.10.07 V.Ivanchenko reorganized methods, add cut dependence in scattering off e- 
// 09.06.08 V.Ivanchenko SelectIsotope is moved to the base class
// 16.06.09 Consolandi rows 109, 111-112, 183, 185-186
//
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
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4CoulombScatteringModel::G4CoulombScatteringModel(const G4String& nam)
  : G4eCoulombScatteringModel(nam)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScatteringModel::~G4CoulombScatteringModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CoulombScatteringModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy, 
				G4double Z, 
				G4double, 
				G4double cutEnergy,
				G4double)
{
  //G4cout << "### G4CoulombScatteringModel::ComputeCrossSectionPerAtom  for " 
  //	 << p->GetParticleName()<<" Z= "<<Z<<" e(MeV)= "<< kinEnergy/MeV 
  //	 <<" cut(MeV)= " << cutEnergy<< G4endl; 
  SetupParticle(p);
  if(kinEnergy < lowEnergyLimit) return 0.0;
  SetupKinematic(kinEnergy, cutEnergy);

  // save lab system kinematics
  G4double xtkin = tkin;
  G4double xmom2 = mom2;
  G4double xinvb = invbeta2;

  // CM system
  iz            = G4int(Z);
  G4double m2   = fNistManager->GetAtomicMassAmu(iz)*amu_c2;
  G4double etot = tkin + mass;
  G4double ptot = sqrt(mom2);

  G4double m12  = mass*mass;

  // 03.09.2009 C.Consaldi suggested to use relativistic reduced mass
  //            from publucation
  // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179
  G4double Ecm=sqrt(m12 + m2*m2 + 2.0*etot*m2);
  G4double mu_rel=mass*m2/Ecm;

  G4double momCM= ptot*m2/Ecm;
  mom2 = momCM*momCM;

  invbeta2 = 1.0 +  mu_rel*mu_rel/mom2; 
  tkin = momCM*sqrt(invbeta2) - mu_rel;
  //

  SetupTarget(Z, tkin);

  G4double xsec = CrossSectionPerAtom();

  // restore Lab system kinematics
  tkin = xtkin;
  mom2 = xmom2;
  invbeta2 = xinvb;
  /*
  G4cout << "e(MeV)= " << kinEnergy/MeV << " xsec(b)= " << xsec/barn  
	 << "cosTetMinNuc= " << cosTetMinNuc
	 << " cosTetMaxNuc= " << cosTetMaxNuc
	 << " cosTetMaxElec= " << cosTetMaxElec
	 << " screenZ= " << screenZ
	 << " formfactA= " << formfactA << G4endl;
  */
  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScatteringModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double cutEnergy, 
			       G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  if(kinEnergy < lowEnergyLimit) return;
  DefineMaterial(couple);
  SetupParticle(dp->GetDefinition());
  SetupKinematic(kinEnergy, cutEnergy);

  // Choose nucleus
  currentElement = SelectRandomAtom(couple,particle,
				    kinEnergy,ecut,kinEnergy);

  G4double Z  = currentElement->GetZ();
  iz          = G4int(Z);
  G4int ia    = SelectIsotopeNumber(currentElement);
  G4double m2 = theParticleTable->GetIonTable()->GetNucleusMass(iz, ia);

  // CM system
  G4double etot = tkin + mass;
  G4double ptot = sqrt(mom2);

  G4double m12  = mass*mass;

  // 03.09.2009 C.Consaldi suggested to use relativistic reduced mass
  //            from publucation
  // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179
  G4double Ecm=sqrt(m12 + m2*m2 + 2.0*etot*m2);
  G4double mu_rel=mass*m2/Ecm;

  G4double momCM= ptot*m2/Ecm;
  mom2 = momCM*momCM;

  invbeta2 = 1.0 +  mu_rel*mu_rel/mom2; 
  tkin = momCM*sqrt(invbeta2) - mu_rel;
  //

  // sample scattering angle in CM system
  SetupTarget(Z, tkin);

  G4double z1 = SampleCosineTheta();
  if(z1 <= 0.0) { return; }
  G4double cost = 1.0 - z1;

  G4double sint = sqrt(z1*(1.0 + cost));
  G4double phi  = twopi * G4UniformRand();

  // kinematics in the Lab system
  G4double bet  = ptot/(etot + m2);
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double pzCM = momCM*cost;
  G4double eCM  = sqrt(mom2 + m12);

  G4ThreeVector v1(momCM*cos(phi)*sint,momCM*sin(phi)*sint,gam*(pzCM + bet*eCM));
  G4ThreeVector dir = dp->GetMomentumDirection(); 
  G4ThreeVector newDirection = v1.unit();
  newDirection.rotateUz(dir);   
  fParticleChange->ProposeMomentumDirection(newDirection);   

  G4double elab = gam*(eCM + bet*pzCM);
 
  // recoil
  G4double finalT = elab - mass;
  G4double trec = kinEnergy - finalT;

  if(finalT <= lowEnergyLimit) { 
    trec = kinEnergy;  
    finalT = 0.0;
  } 
    
  fParticleChange->SetProposedKineticEnergy(finalT);

  //  G4cout << "sint= " << sint << " Erec(eV)= " << erec/eV << G4endl;

  G4double tcut = recoilThreshold;
  if(pCuts) { tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]); } 
  /*  
  G4cout << "sint= " << sint << " Erec(eV)= " << erec/eV
	 << " tcut(eV)= " << tcut/eV << " th(eV)= " << recoilThreshold/eV
	 << " cut(eV)= " << (*pCuts)[currentMaterialIndex]/eV
	 << "  "  << fvect->size()
	 << G4endl;
  */
  if(trec > tcut) {
    G4ParticleDefinition* ion = theParticleTable->FindIon(iz, ia, 0, iz);
    G4double plab = sqrt(finalT*(finalT + 2.0*mass));
    G4ThreeVector p2 = (ptot*dir - plab*newDirection).unit();
    G4DynamicParticle* newdp  = new G4DynamicParticle(ion, p2, trec);
    fvect->push_back(newdp);
  } else if(trec > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(trec);
    fParticleChange->ProposeNonIonizingEnergyDeposit(trec);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

