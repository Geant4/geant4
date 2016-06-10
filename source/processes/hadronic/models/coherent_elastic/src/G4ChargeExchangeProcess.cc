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
// $Id: G4ChargeExchangeProcess.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
//
// Geant4 Hadron Charge Exchange Process -- source file
//
// Created 21 April 2006 V.Ivanchenko
//
// Modified:
// 24-Apr-06 V.Ivanchenko add neutron scattering on hydrogen from CHIPS
// 07-Jun-06 V.Ivanchenko fix problem of rotation of final state
// 25-Jul-06 V.Ivanchenko add 19 MeV low energy for CHIPS
// 23-Jan-07 V.Ivanchenko add cross section interfaces with Z and A
//                        and do not use CHIPS for cross sections
// 14-Sep-12 M.Kelsey -- Pass subType code to base ctor
// 06-Aug-15 A.Ribon migrating to G4Pow

#include "G4ChargeExchangeProcess.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4PhysicsLinearVector.hh"

#include "G4Pow.hh"


G4ChargeExchangeProcess::G4ChargeExchangeProcess(const G4String& procName)
  : G4HadronicProcess(procName,fChargeExchange), first(true)
{
  thEnergy = 20.*MeV;
  pPDG = 0;
  verboseLevel= 1;
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
  theHe3      = G4He3::He3();
}

G4ChargeExchangeProcess::~G4ChargeExchangeProcess()
{
  if (factors) delete factors;
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
      factors = new G4PhysicsLinearVector(0.0,2.0*GeV,n);
      for(size_t i=0; i<n; i++) {factors->PutValue(i,F[i]);}

    } else {

      G4double F[n] = {0.50,0.45,0.40,0.35,0.30,0.25,0.06,0.04,0.005,0.0};
      factors = new G4PhysicsLinearVector(0.0,4.0*GeV,n);
      for(size_t i=0; i<n; i++) {factors->PutValue(i,F[i]);}
    }
    //factors->SetSpline(true);

    if(verboseLevel>1)
      G4cout << "G4ChargeExchangeProcess for "
	     << theParticle->GetParticleName()
	     << G4endl;
  }
  G4HadronicProcess::BuildPhysicsTable(aParticleType);
}

G4double G4ChargeExchangeProcess::GetElementCrossSection(
                                  const G4DynamicParticle* dp,
				  const G4Element* elm,
				  const G4Material* mat)
{
  // gives the microscopic cross section in GEANT4 internal units
  G4double Z = elm->GetZ();
  G4int iz = G4int(Z);
  G4double x = 0.0;

  // The process is effective only above the threshold
  if(iz == 1 || dp->GetKineticEnergy() < thEnergy) return x;

  if(verboseLevel>1)
    G4cout << "G4ChargeExchangeProcess compute GHAD CS for element "
	   << elm->GetName()
	   << G4endl;
  x = store->GetCrossSection(dp, elm, mat);

  if(verboseLevel>1)
    G4cout << "G4ChargeExchangeProcess cross(mb)= " << x/millibarn
           << "  E(MeV)= " << dp->GetKineticEnergy()
	   << "  " << theParticle->GetParticleName()
           << "  in Z= " << iz
	   << G4endl;
  G4bool b;
  G4double A = elm->GetN();
  G4double ptot = dp->GetTotalMomentum();
  x *= factors->GetValue(ptot, b)/G4Pow::GetInstance()->powA(A, 0.42);
  if(theParticle == thePiPlus || theParticle == theProton ||
     theParticle == theKPlus  || theParticle == theANeutron)
    { x *= (1.0 - Z/A); }

  else if(theParticle == thePiMinus || theParticle == theNeutron ||
          theParticle == theKMinus  || theParticle == theAProton)
    { x *= Z/A; }

  if(theParticle->GetPDGMass() < GeV) {
    if(ptot > 2.*GeV) x *= 4.0*GeV*GeV/(ptot*ptot);
  }

  if(verboseLevel>1) 
    G4cout << "Corrected cross(mb)= " << x/millibarn << G4endl;

  return x;
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
