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
//
// $Id: G4ChargeExchangeProcess.cc,v 1.9 2007/01/30 10:23:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-03 $
//
//
// Geant4 Hadron Elastic Scattering Process -- header file
//
// Created 21 April 2006 V.Ivanchenko
//
// Modified:
// 24-Apr-06 V.Ivanchenko add neutron scattering on hydrogen from CHIPS
// 07-Jun-06 V.Ivanchenko fix problem of rotation of final state
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy for CHIPS
// 23-Jan-07 V.Ivanchenko add cross section interfaces with Z and A
//                        and do not use CHIPS for cross sections
//

#include "G4ChargeExchangeProcess.hh"
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
#include "G4HadronElastic.hh"
#include "G4PhysicsLinearVector.hh"

G4ChargeExchangeProcess::G4ChargeExchangeProcess(const G4String& procName)
  : G4HadronicProcess(procName), first(true)
{
  thEnergy = 19.*MeV;
  verboseLevel= 1;
  qCManager = 0;
  AddDataSet(new G4HadronElasticDataSet);
  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theAProton  = G4AntiProton::AntiProton();
  theANeutron = G4AntiNeutron::AntiNeutron();
  thePiPlus   = G4PionPlus::PionPlus();
  thePiMinus  = G4PionMinus::PionMinus();
  thePiZero   = G4PionZero::PionZero();
  theKPlus    = G4KaonPlus::KaonPlus();
  theKMinus   = G4KaonMinus::KaonMinus();
  theK0S      = G4KaonZeroShort::KaonZeroShort();
  theK0L      = G4KaonZeroLong::KaonZeroLong();
  theL        = G4Lambda::Lambda();
  theAntiL    = G4AntiLambda::AntiLambda();
  theSPlus    = G4SigmaPlus::SigmaPlus();
  theASPlus   = G4AntiSigmaPlus::AntiSigmaPlus();
  theSMinus   = G4SigmaMinus::SigmaMinus();
  theASMinus  = G4AntiSigmaMinus::AntiSigmaMinus();
  theS0       = G4SigmaZero::SigmaZero();
  theAS0      = G4AntiSigmaZero::AntiSigmaZero();
  theXiMinus  = G4XiMinus::XiMinus();
  theXi0      = G4XiZero::XiZero();
  theAXiMinus = G4AntiXiMinus::AntiXiMinus();
  theAXi0     = G4AntiXiZero::AntiXiZero();
  theOmega    = G4OmegaMinus::OmegaMinus();
  theAOmega   = G4AntiOmegaMinus::AntiOmegaMinus();
  theD        = G4Deuteron::Deuteron();
  theT        = G4Triton::Triton();
  theA        = G4Alpha::Alpha();
  theA        = G4He3::He3();
}

G4ChargeExchangeProcess::~G4ChargeExchangeProcess()
{
  delete factors;
}

void G4ChargeExchangeProcess::SetQElasticCrossSection(G4VQCrossSection* p)
{
  qCManager = p;
}

void G4ChargeExchangeProcess::
BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(first) {
    first = false;
    theParticle = &aParticleType;
    pPDG = theParticle->GetPDGEncoding();

    store = G4HadronicProcess::GetCrossSectionDataStore();

    const size_t n = 10;
    if(theParticle == thePiPlus || theParticle == thePiMinus ||
       theParticle == theKPlus  || theParticle == theKMinus ||
       theParticle == theK0S    || theParticle == theK0L) {

      G4double F[n] = {0.33,0.27,0.29,0.31,0.27,0.18,0.13,0.1,0.09,0.07};
      factors = new G4PhysicsLinearVector(0.0,1.8*GeV,n);
      for(size_t i=0; i<n; i++) {factors->PutValue(i,F[i]);}

    } else {

      G4double F[n] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
      factors = new G4PhysicsLinearVector(0.0,3.6*GeV,n);
      for(size_t i=0; i<n; i++) {factors->PutValue(i,F[i]);}
    }

    if(verboseLevel>1)
      G4cout << "G4ChargeExchangeProcess for "
	     << theParticle->GetParticleName()
	     << G4endl;
  }
  store->BuildPhysicsTable(aParticleType);
}

G4double G4ChargeExchangeProcess::GetMeanFreePath(const G4Track& track,
						  G4double,
						  G4ForceCondition* cond)
{
  *cond = NotForced;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  const G4Material* material = track.GetMaterial();
  cross = 0.0;
  G4double x = DBL_MAX;

  // The process is effective only above the threshold
  if(dp->GetKineticEnergy() < thEnergy) return x;

  // Compute cross sesctions
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  G4double temp = material->GetTemperature();
  G4int nelm   = material->GetNumberOfElements();
  if(verboseLevel>1)
    G4cout << "G4ChargeExchangeProcess get mfp for "
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
    G4cout << "G4ChargeExchangeProcess cross(1/mm)= " << cross
           << "  E(MeV)= " << dp->GetKineticEnergy()
	   << "  " << theParticle->GetParticleName()
           << "  in " << material->GetName()
	   << G4endl;
  if(cross > DBL_MIN) x = 1./cross;

  return x;
}

G4double G4ChargeExchangeProcess::GetMicroscopicCrossSection(
                                  const G4DynamicParticle* dp,
				  const G4Element* elm,
				  G4double temp)
{
  // gives the microscopic cross section in GEANT4 internal units
  G4double Z = elm->GetZ();
  G4int iz = G4int(Z);
  G4double x = 0.0;
  if(iz == 1) return x;

  if(verboseLevel>1)
    G4cout << "G4ChargeExchangeProcess compute GHAD CS for element "
	   << elm->GetName()
	   << G4endl;
  x = store->GetCrossSection(dp, elm, temp);

  // NaN finder
  if(!(x < 0.0 || x >= 0.0)) {
    if (verboseLevel > -1) {
      G4cout << "G4ChargeExchangeProcess WARNING: Z= " << iz  
	     << " pdg= " <<  pPDG
	     << " mom(GeV)= " << dp->GetTotalMomentum()/GeV 
	     << " cross= " << x 
	     << " set to zero"
	     << G4endl; 
    }
    x = 0.0;
  }

  if(verboseLevel>1)
    G4cout << "G4ChargeExchangeProcess cross(mb)= " << x/millibarn
           << "  E(MeV)= " << dp->GetKineticEnergy()
	   << "  " << theParticle->GetParticleName()
           << "  in Z= " << iz
	   << G4endl;
  G4bool b;
  G4double A = elm->GetN();
  x *= factors->GetValue(dp->GetTotalMomentum(), b)/std::pow(A, 0.42);
  if(theParticle == thePiPlus || theParticle == theProton ||
     theParticle == theKPlus  || theParticle == theANeutron)
    x *= (1.0 - Z/A);

  else if(theParticle == thePiMinus || theParticle == theNeutron ||
          theParticle == theKMinus  || theParticle == theAProton)
    x *= Z/A;

  return x;
}

G4VParticleChange* G4ChargeExchangeProcess::PostStepDoIt(
				  const G4Track& track,
				  const G4Step& step)
{
  G4ForceCondition* cn = 0;
  aParticleChange.Initialize(track);
  G4double mfp = GetMeanFreePath(track, 0.0, cn);
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
  G4double A = G4double(G4int(elm->GetN()+0.5));

  // Select isotope
  G4IsotopeVector* isv = elm->GetIsotopeVector(); 
  G4int ni = 0;
  if(isv) ni = isv->size();

  if(ni == 1) { 
    A = G4double((*isv)[0]->GetN());
  } else if(ni > 1) {

    G4double* ab = elm->GetRelativeAbundanceVector();
    G4int j = -1;
    ni--;
    G4double y = G4UniformRand();
    do {
      j++;
      y -= ab[j];
    } while (y > 0.0 && j < ni);
    A = G4double((*isv)[j]->GetN());
  } 
  G4HadronicInteraction* hadi =
    ChooseHadronicInteraction( kineticEnergy, material, elm);

  // Initialize the hadronic projectile from the track
  G4HadProjectile thePro(track);
  if(verboseLevel>1)
    G4cout << "G4ChargeExchangeProcess::PostStepDoIt for "
	   << theParticle->GetParticleName()
	   << " Target Z= " << Z
	   << " A= " << A << G4endl;
  targetNucleus.SetParameters(A, Z);

  aParticleChange.Initialize(track);
  G4HadFinalState* result = hadi->ApplyYourself(thePro, targetNucleus);
  G4ThreeVector indir = track.GetMomentumDirection();
  G4int nsec = result->GetNumberOfSecondaries();

  if(verboseLevel>1)
    G4cout << "Efin= " << result->GetEnergyChange()
	   << " de= " << result->GetLocalEnergyDeposit()
	   << " nsec= " << nsec
	   << G4endl;


  if(nsec > 0) {
    aParticleChange.ProposeEnergy(0.0);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.ProposeLocalEnergyDeposit(result->GetLocalEnergyDeposit());
    aParticleChange.SetNumberOfSecondaries(nsec);
    for(G4int j=0; j<nsec; j++) {
      G4DynamicParticle* p = result->GetSecondary(j)->GetParticle();
      G4ThreeVector pdir = p->GetMomentumDirection();
      // G4cout << "recoil " << pdir << G4endl;
      pdir = pdir.rotateUz(indir);
      // G4cout << "recoil rotated " << pdir << G4endl;
      p->SetMomentumDirection(pdir);
      aParticleChange.AddSecondary(p);
    }
  }
  result->Clear();

  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

G4bool G4ChargeExchangeProcess::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
  const G4ParticleDefinition* p = &aParticleType;
  return (p == thePiPlus || p == thePiMinus ||
          p == theProton || p == theNeutron ||
          p == theAProton|| p == theANeutron||
	  p == theKPlus  || p == theKMinus  ||
	  p == theK0S    || p == theK0L     ||
	  p == theL);
}

void G4ChargeExchangeProcess::
DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  store->DumpPhysicsTable(aParticleType);
}
