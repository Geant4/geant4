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
// $Id: G4UHadronElasticProcess.cc,v 1.18 2006/07/12 17:15:41 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
// Geant4 Hadron Elastic Scattering Process -- header file
// 
// Created 21 April 2006 V.Ivanchenko
//  
// Modified:
// 24-Apr-06 V.Ivanchenko add neutron scattering on hydrogen from CHIPS
// 07-Jun-06 V.Ivanchenko fix problem of rotation of final state
//
//

#include "G4UHadronElasticProcess.hh"
#include "globals.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4VQCrossSection.hh"
#include "G4QElasticCrossSection.hh"
#include "G4QCHIPSWorld.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4HadronElastic.hh"
 
G4UHadronElasticProcess::G4UHadronElasticProcess(const G4String& pName, G4bool fl)
  : G4HadronicProcess(pName), flagHP(fl), first(true)
{
  AddDataSet(new G4HadronElasticDataSet);
  theProton = G4Proton::Proton();
  theNeutron = G4Neutron::Neutron();
  thEnergy = 1.*keV;
  verboseLevel= 1;
  qCManager = 0;
}

G4UHadronElasticProcess::~G4UHadronElasticProcess()
{
}

void G4UHadronElasticProcess::SetQElasticCrossSection(G4VQCrossSection* p)
{
  qCManager = p;
}

void G4UHadronElasticProcess::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(first) {
    first = false;
    if(!qCManager) qCManager = G4QElasticCrossSection::GetPointer();
    theParticle = &aParticleType;
    pPDG = theParticle->GetPDGEncoding();
    if(theParticle == theNeutron && flagHP) 
      AddDataSet(new G4NeutronHPElasticData());

    store = G4HadronicProcess::GetCrossSectionDataStore();
     
    if(verboseLevel>1) 
      G4cout << "G4UHadronElasticProcess for " 
	     << theParticle->GetParticleName() 
	     << G4endl; 
  }
  store->BuildPhysicsTable(aParticleType);
}

G4double G4UHadronElasticProcess::GetMeanFreePath(const G4Track& track, 
						  G4double, 
						  G4ForceCondition* cond)
{
  *cond = NotForced;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  const G4Material* material = track.GetMaterial();
  cross = 0.0;
  G4double x = DBL_MAX;

  // The process is effective only above the threshold
  //  if(dp->GetKineticEnergy() < thEnergy) return x;

  // Compute cross sesctions
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  G4double temp = material->GetTemperature();
  G4int nelm   = material->GetNumberOfElements();
  xsecH[0] = 0.0;
  xsecH[1] = 0.0;
  if(verboseLevel>1) 
    G4cout << "G4UHadronElasticProcess get mfp for " 
	   << theParticle->GetParticleName() 
	   << "  p(GeV)= " << dp->GetTotalMomentum()/GeV
	   << " in " << material->GetName()
	   << G4endl; 
 
  for (G4int i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double x = GetMicroscopicCrossSection(dp, elm, temp);
    cross += theAtomNumDensityVector[i]*x;
    xsec[i] = cross;
  }
  if(verboseLevel>1) 
    G4cout << "G4UHadronElasticProcess cross(1/mm)= " << cross 
           << "  E(MeV)= " << dp->GetKineticEnergy()
	   << "  " << theParticle->GetParticleName()
           << "  in " << material->GetName()
	   << G4endl;
  if(cross > DBL_MIN) x = 1./cross;

  return x;
}

G4double G4UHadronElasticProcess::GetMicroscopicCrossSection(
                                  const G4DynamicParticle* dp, 
				  const G4Element* elm, 
				  G4double temp)
{
  // gives the microscopic cross section in GEANT4 internal units
  G4int iz = G4int(elm->GetZ());
  G4double x = 0.0;
  // CHIPS cross sections
  if(iz <= 2 && (theParticle == theProton || theParticle == theNeutron)) {
    G4double momentum = dp->GetTotalMomentum();
    if(iz == 1) {
      G4IsotopeVector* isv = elm->GetIsotopeVector(); 
      G4int ni = 0;
      if(isv) ni = isv->size();
      if(ni > 0) {
	G4double* ab = elm->GetRelativeAbundanceVector();
	x = 0.0;
	for(G4int j=0; j<ni; j++) {
	  G4int N = elm->GetIsotope(j)->GetN() - 1;
	  if(N == 0 || N == 1) {
	    if(verboseLevel>1) 
	      G4cout << "G4UHadronElasticProcess compute CHIPS CS for Z= 1, N= " 
		     << N << " pdg= " << pPDG 
		     << " mom(GeV)= " << momentum/GeV 
		     << "  " << qCManager << G4endl; 
	    G4double y = ab[j]*
	      qCManager->GetCrossSection(false,momentum,1,N,pPDG);
	    xsecH[N] += y;
	    x += y;
	  }
	}
      } else {
	if(verboseLevel>1) 
	  G4cout << "G4UHadronElasticProcess compute CHIPS CS for Z= 1, N= 0" 
		 << " pdg= " << pPDG 
		 << " mom(GeV)= " << momentum/GeV << "  " << qCManager << G4endl; 
	x = qCManager->GetCrossSection(false,momentum,1,0,pPDG);
	xsecH[0] = x;
      }
    } else {
      if(verboseLevel>1) 
	G4cout << "G4UHadronElasticProcess compute CHIPS CS for Z= 2, N=2 " 
	       << G4endl; 
      x = qCManager->GetCrossSection(false,momentum,2,2,pPDG);
    }
  } else {
    if(verboseLevel>1) 
      G4cout << "G4UHadronElasticProcess compute GHAD CS for element " 
	     << elm->GetName() 
	     << G4endl; 
    x = store->GetCrossSection(dp, elm, temp);
  }
  if(verboseLevel>1) 
    G4cout << "G4UHadronElasticProcess cross(mb)= " << x/millibarn 
           << "  E(MeV)= " << dp->GetKineticEnergy()
	   << "  " << theParticle->GetParticleName()
           << "  in Z= " << iz
	   << G4endl;

  return x;
}

G4VParticleChange* G4UHadronElasticProcess::PostStepDoIt(
				  const G4Track& track, 
				  const G4Step& step)
{
  G4ForceCondition   cn;
  aParticleChange.Initialize(track);
  G4double mfp = GetMeanFreePath(track, 0.0, &cn);
  if(mfp == DBL_MAX) return G4VDiscreteProcess::PostStepDoIt(track,step);

  G4double kineticEnergy = track.GetKineticEnergy();
  G4Material* material = track.GetMaterial();

  // Select element
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4Element* elm = (*theElementVector)[0];
  G4int nelm = material->GetNumberOfElements() - 1;
  if (nelm > 0) {
    G4double x = G4UniformRand()*cross;
    G4int i = -1;
    do {i++;} while (x > xsec[i] && i < nelm);
    elm = (*theElementVector)[i];
  }
  G4double Z = elm->GetZ();
  G4double A = elm->GetN();
  G4int iz = G4int(Z);

  // Select isotope
  G4IsotopeVector* isv = elm->GetIsotopeVector(); 
  G4int ni = 0;
  if(isv) ni = isv->size();
  if(ni == 1) { 
    A = G4double(elm->GetIsotope(0)->GetN());
  } else if(ni > 1) {
    if(iz == 1 && theParticle == theProton) {
      A = 1.;
      if(G4UniformRand()*(xsecH[0] + xsec[1]) > xsec[0]) A = 2.;
    } else {
      G4double* ab = elm->GetRelativeAbundanceVector();
      G4double y = G4UniformRand();
      G4int j = -1;
      ni--;
      do {
	j++;
	y -= ab[j];
      } while (y > 0.0 && j < ni);
      A = G4double(elm->GetIsotope(j)->GetN());
    }
  }
  G4HadronicInteraction* hadi = 
    ChooseHadronicInteraction( kineticEnergy, material, elm);

  // Initialize the hadronic projectile from the track
  //  G4cout << "track " << track.GetDynamicParticle()->Get4Momentum()<<G4endl;
  G4HadProjectile thePro(track);
  if(verboseLevel>1) 
    G4cout << "G4UHadronElasticProcess::PostStepDoIt for " 
	   << theParticle->GetParticleName() 
	   << " Target Z= " << Z 
	   << " A= " << A << G4endl; 
  targetNucleus.SetParameters(A, Z);

  aParticleChange.Initialize(track);
  G4HadFinalState* result = hadi->ApplyYourself(thePro, targetNucleus);
  G4ThreeVector indir = track.GetMomentumDirection();
  G4ThreeVector outdir = (result->GetMomentumChange()).rotateUz(indir);
  
  if(verboseLevel>1) 
    G4cout << "Efin= " << result->GetEnergyChange()
	   << " de= " << result->GetLocalEnergyDeposit()
	   << " nsec= " << result->GetNumberOfSecondaries()
	   << " dir= " << outdir
	   << G4endl;
  
  aParticleChange.ProposeEnergy(result->GetEnergyChange());
  aParticleChange.ProposeMomentumDirection(outdir);
  if(result->GetNumberOfSecondaries() > 0) {
    aParticleChange.SetNumberOfSecondaries(1);
    G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();
    G4ThreeVector pdir = p->GetMomentumDirection();
    // G4cout << "recoil " << pdir << G4endl;
    pdir = pdir.rotateUz(indir);
    // G4cout << "recoil rotated " << pdir << G4endl;
    p->SetMomentumDirection(pdir);
    aParticleChange.AddSecondary(p);
  } else {
    aParticleChange.SetNumberOfSecondaries(0);
    aParticleChange.ProposeLocalEnergyDeposit(result->GetLocalEnergyDeposit());
  }
  result->Clear();

  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

G4bool G4UHadronElasticProcess::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return (aParticleType == *(G4PionPlus::PionPlus()) ||
           aParticleType == *(G4PionMinus::PionMinus()) ||
           aParticleType == *(G4KaonPlus::KaonPlus()) ||
           aParticleType == *(G4KaonZeroShort::KaonZeroShort()) ||
           aParticleType == *(G4KaonZeroLong::KaonZeroLong()) ||
           aParticleType == *(G4KaonMinus::KaonMinus()) ||
           aParticleType == *(G4Proton::Proton()) ||
           aParticleType == *(G4AntiProton::AntiProton()) ||
           aParticleType == *(G4Neutron::Neutron()) ||
           aParticleType == *(G4AntiNeutron::AntiNeutron()) ||
           aParticleType == *(G4Lambda::Lambda()) ||
           aParticleType == *(G4AntiLambda::AntiLambda()) ||
           aParticleType == *(G4SigmaPlus::SigmaPlus()) ||
           aParticleType == *(G4SigmaZero::SigmaZero()) ||
           aParticleType == *(G4SigmaMinus::SigmaMinus()) ||
           aParticleType == *(G4AntiSigmaPlus::AntiSigmaPlus()) ||
           aParticleType == *(G4AntiSigmaZero::AntiSigmaZero()) ||
           aParticleType == *(G4AntiSigmaMinus::AntiSigmaMinus()) ||
           aParticleType == *(G4XiZero::XiZero()) ||
           aParticleType == *(G4XiMinus::XiMinus()) ||
           aParticleType == *(G4AntiXiZero::AntiXiZero()) ||
           aParticleType == *(G4AntiXiMinus::AntiXiMinus()) ||
           aParticleType == *(G4Deuteron::Deuteron()) ||
           aParticleType == *(G4Triton::Triton()) ||
           aParticleType == *(G4Alpha::Alpha()) ||
           aParticleType == *(G4OmegaMinus::OmegaMinus()) ||
           aParticleType == *(G4AntiOmegaMinus::AntiOmegaMinus()));
}

void G4UHadronElasticProcess::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  store->DumpPhysicsTable(aParticleType);
}
