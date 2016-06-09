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
// $Id: G4CoulombScatteringModel.cc,v 1.49 2010-05-27 14:22:05 vnivanch Exp $
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
//
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.10.06 V.Ivanchenko use inheritance from G4eCoulombScatteringModel
// 09.10.07 V.Ivanchenko reorganized methods, add cut dependence in scattering off e- 
// 09.06.08 V.Ivanchenko SelectIsotope is moved to the base class
// 16.06.09 Consolandi rows 109, 111-112, 183, 185-186
// 27.05.10 V.Ivanchenko added G4WentzelOKandVIxSection class to
//              compute cross sections and sample scattering angle
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
#include "G4ProcessManager.hh"
#include "G4NucleiProperties.hh"

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
  G4double xsec = 0.0;
  if(p != particle) { SetupParticle(p); }
  if(kinEnergy <= 0.0) { return 0.0; }
  DefineMaterial(CurrentCouple());

  // Lab system
  G4int iz = G4int(Z);
  G4double etot = kinEnergy + mass;
  G4double m2 = fNistManager->GetAtomicMassAmu(iz)*amu_c2;

  // 03.09.2009 C.Consaldi suggested to use relativistic reduced mass
  //            from publucation
  // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179
  G4double Ecm  = sqrt(mass*mass + m2*m2 + 2.0*etot*m2);
  G4double mu_rel = mass*m2/Ecm;
  G4double tkin = Ecm - mu_rel;
  wokvi->SetRelativisticMass(mu_rel);
  
  cosTetMinNuc = wokvi->SetupKinematic(tkin, currentMaterial);
  if(cosThetaMax < cosTetMinNuc) {
    cosTetMinNuc = wokvi->SetupTarget(iz, cutEnergy);
    cosTetMaxNuc = cosThetaMax; 
    if(iz == 1 && cosTetMaxNuc < 0.0 && particle == theProton) { 
      cosTetMaxNuc = 0.0; 
    }
    xsec =  wokvi->ComputeNuclearCrossSection(cosTetMinNuc, cosTetMaxNuc);
    elecRatio = wokvi->ComputeElectronCrossSection(cosTetMinNuc, cosThetaMax);
    xsec += elecRatio;
    if(xsec > 0.0) { elecRatio /= xsec; }  
  }
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
  //  if(kinEnergy < lowEnergyLimit) { return; }
  if(kinEnergy < lowEnergyLimit) {
    fParticleChange->SetProposedKineticEnergy(0.0);
    fParticleChange->ProposeLocalEnergyDeposit(kinEnergy);
    fParticleChange->ProposeNonIonizingEnergyDeposit(kinEnergy);
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange->ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange->ProposeTrackStatus(fStopAndKill); }
    return;
  }
  DefineMaterial(couple);
  SetupParticle(dp->GetDefinition());

  // Choose nucleus
  currentElement = SelectRandomAtom(couple,particle,
				    kinEnergy,cutEnergy,kinEnergy);

  G4double Z = currentElement->GetZ();
  G4int iz = G4int(Z);
  G4int ia = SelectIsotopeNumber(currentElement);
  G4double targetMass = G4NucleiProperties::GetNuclearMass(ia, iz);
  
  if(ComputeCrossSectionPerAtom(particle,kinEnergy, Z,
				kinEnergy, cutEnergy, kinEnergy) == 0.0) 
    { return; }

  G4ThreeVector newDirection = 
    wokvi->SampleSingleScattering(cosTetMinNuc, cosTetMaxNuc, elecRatio);

  // kinematics in the Lab system
  G4double etot = mass + kinEnergy;
  G4double ptot = sqrt(kinEnergy*(etot + mass));
  G4double bet  = ptot/(etot + targetMass);
  G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));
  G4double eCM  = sqrt(mass*mass + targetMass*targetMass + 2*targetMass*etot);
  G4double pCM  = ptot*targetMass/eCM;
  G4double e1   = sqrt(mass*mass + pCM*pCM);

  newDirection *= pCM;

  G4ThreeVector v1(newDirection.x(),newDirection.y(),gam*(newDirection.z() + bet*e1));
  G4double finalT = gam*(e1 + bet*newDirection.z()) - mass;
  newDirection = v1.unit();

  G4ThreeVector dir = dp->GetMomentumDirection(); 
  newDirection.rotateUz(dir);   
  fParticleChange->ProposeMomentumDirection(newDirection);   
 
  // recoil
  G4double trec = kinEnergy - finalT;
  if(finalT <= lowEnergyLimit) { 
    trec = kinEnergy;  
    finalT = 0.0;
  } 
    
  if(finalT <= lowEnergyLimit) { 
    trec = kinEnergy;  
    finalT = 0.0;
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange->ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange->ProposeTrackStatus(fStopAndKill); }
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

