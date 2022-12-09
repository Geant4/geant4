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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmExtraParameters
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2019
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmExtraParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4EmExtraParametersMessenger.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmExtraParameters::G4EmExtraParameters()
{
  theMessenger = new G4EmExtraParametersMessenger(this);
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmExtraParameters::~G4EmExtraParameters()
{
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmExtraParameters::Initialise()
{
  quantumEntanglement = false;
  directionalSplitting = false;
  directionalSplittingTarget.set(0.,0.,0.);
  directionalSplittingRadius = 0.;

  dRoverRange = 0.2;
  finalRange = CLHEP::mm;
  dRoverRangeMuHad = 0.2;
  finalRangeMuHad = 0.1*CLHEP::mm;
  dRoverRangeLIons = 0.2;
  finalRangeLIons = 0.1*CLHEP::mm;
  dRoverRangeIons = 0.2;
  finalRangeIons = 0.1*CLHEP::mm;

  m_regnamesForced.clear();
  m_procForced.clear();
  m_lengthForced.clear();
  m_weightForced.clear();
  m_regnamesSubCut.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....


void G4EmExtraParameters::PrintWarning(G4ExceptionDescription& ed) const
{
  G4Exception("G4EmExtraParameters", "em0044", JustWarning, ed);
}

G4String G4EmExtraParameters::CheckRegion(const G4String& reg) const
{
  G4String r = reg;
  if(r == "" || r == "world" || r == "World") {
    r = "DefaultRegionForTheWorld";
  }
  return r;
}

void G4EmExtraParameters::SetStepFunction(G4double v1, G4double v2)
{
  if(v1 > 0.0 && v1 <= 1.0 && v2 > 0.0) {
    dRoverRange = v1;
    finalRange = v2;
  } else {
    G4ExceptionDescription ed;
    ed << "Values of step function are out of range: " 
       << v1 << ", " << v2/CLHEP::mm << " mm - are ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmExtraParameters::GetStepFunctionP1() const
{
  return dRoverRange;
}

G4double G4EmExtraParameters::GetStepFunctionP2() const
{
  return finalRange;
}

void G4EmExtraParameters::SetStepFunctionMuHad(G4double v1, G4double v2)
{
  if(v1 > 0.0 && v1 <= 1.0 && v2 > 0.0) {
    dRoverRangeMuHad = v1;
    finalRangeMuHad = v2;
  } else {
    G4ExceptionDescription ed;
    ed << "Values of step function are out of range: " 
       << v1 << ", " << v2/CLHEP::mm << " mm - are ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmExtraParameters::GetStepFunctionMuHadP1() const
{
  return dRoverRangeMuHad;
}

G4double G4EmExtraParameters::GetStepFunctionMuHadP2() const
{
  return finalRangeMuHad;
}

void G4EmExtraParameters::SetStepFunctionLightIons(G4double v1, G4double v2)
{
  if(v1 > 0.0 && v1 <= 1.0 && v2 > 0.0) {
    dRoverRangeLIons = v1;
    finalRangeLIons = v2;
  } else {
    G4ExceptionDescription ed;
    ed << "Values of step function are out of range: " 
       << v1 << ", " << v2/CLHEP::mm << " mm - are ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmExtraParameters::GetStepFunctionLightIonsP1() const
{
  return dRoverRangeLIons;
}

G4double G4EmExtraParameters::GetStepFunctionLightIonsP2() const
{
  return finalRangeLIons;
}

void G4EmExtraParameters::SetStepFunctionIons(G4double v1, G4double v2)
{
  if(v1 > 0.0 && v1 <= 1.0 && v2 > 0.0) {
    dRoverRangeIons = v1;
    finalRangeIons = v2;
  } else {
    G4ExceptionDescription ed;
    ed << "Values of step function are out of range: " 
       << v1 << ", " << v2/CLHEP::mm << " mm - are ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmExtraParameters::GetStepFunctionIonsP1() const
{
  return dRoverRangeIons;
}

G4double G4EmExtraParameters::GetStepFunctionIonsP2() const
{
  return finalRangeIons;
}

void G4EmExtraParameters::FillStepFunction(const G4ParticleDefinition* part, G4VEnergyLossProcess* proc) const
{
  // electron and positron
  if (11 == std::abs(part->GetPDGEncoding())) {
    proc->SetStepFunction(dRoverRange, finalRange);

    // all heavy ions
  } else if ("GenericIon" == part->GetParticleName()) {
    proc->SetStepFunction(dRoverRangeIons, finalRangeIons);

    // light nucleus and anti-nucleus
  } else if (part->GetParticleType() == "nucleus" || part->GetParticleType() == "anti_nucleus") { 
    proc->SetStepFunction(dRoverRangeLIons, finalRangeLIons);

    // other particles
  } else {
    proc->SetStepFunction(dRoverRangeMuHad, finalRangeMuHad);
  }
}

void G4EmExtraParameters::AddPAIModel(const G4String& particle,
                                      const G4String& region,
                                      const G4String& type)
{
  G4String r = CheckRegion(region);
  std::size_t nreg =  m_regnamesPAI.size();
  for(std::size_t i=0; i<nreg; ++i) {
    if((m_particlesPAI[i] == particle || 
        m_particlesPAI[i] == "all" || 
        particle == "all") && 
       (m_regnamesPAI[i] == r || 
        m_regnamesPAI[i] == "DefaultRegionForTheWorld" || 
        r == "DefaultRegionForTheWorld") ) {

      m_typesPAI[i] = type;
      if(particle == "all") { m_particlesPAI[i] = particle; }
      if(r == "DefaultRegionForTheWorld") { m_regnamesPAI[i] = r; }
      return;
    }
  }
  m_particlesPAI.push_back(particle);
  m_regnamesPAI.push_back(r);
  m_typesPAI.push_back(type);
}

const std::vector<G4String>& G4EmExtraParameters::ParticlesPAI() const
{
  return m_particlesPAI;
}

const std::vector<G4String>& G4EmExtraParameters::RegionsPAI() const
{
  return m_regnamesPAI;
}

const std::vector<G4String>& G4EmExtraParameters::TypesPAI() const
{
  return m_typesPAI;
}

void G4EmExtraParameters::AddPhysics(const G4String& region, 
                                     const G4String& type)
{
  G4String r = CheckRegion(region);
  std::size_t nreg =  m_regnamesPhys.size();
  for(std::size_t i=0; i<nreg; ++i) {
    if(r == m_regnamesPhys[i]) { return; }
  }
  m_regnamesPhys.push_back(r);
  m_typesPhys.push_back(type);
}

const std::vector<G4String>& G4EmExtraParameters::RegionsPhysics() const
{
  return m_regnamesPhys;
}

const std::vector<G4String>& G4EmExtraParameters::TypesPhysics() const
{
  return m_typesPhys;
}

void G4EmExtraParameters::SetSubCutRegion(const G4String& region)
{
  const G4String& r = CheckRegion(region);
  std::size_t nreg =  m_regnamesSubCut.size();
  for(std::size_t i=0; i<nreg; ++i) {
    if(r == m_regnamesSubCut[i]) { 
      return; 
    }
  }
  m_regnamesSubCut.push_back(r);
}

void 
G4EmExtraParameters::SetProcessBiasingFactor(const G4String& procname, 
                                             G4double val, G4bool wflag)
{
  if(val > 0.0) {
    std::size_t n =  m_procBiasedXS.size();
    for(std::size_t i=0; i<n; ++i) {
      if(procname == m_procBiasedXS[i]) { 
	m_factBiasedXS[i] = val;
	m_weightBiasedXS[i]= wflag;
	return; 
      }
    }
    m_procBiasedXS.push_back(procname);
    m_factBiasedXS.push_back(val);
    m_weightBiasedXS.push_back(wflag);
  } else {
    G4ExceptionDescription ed;
    ed << "Process: " << procname << " XS biasing factor " 
       << val << " is negative - ignored"; 
    PrintWarning(ed);
  }
}

void 
G4EmExtraParameters::ActivateForcedInteraction(const G4String& procname, 
                                               const G4String& region,
                                               G4double length, 
                                               G4bool wflag)
{
  const G4String& r = CheckRegion(region);
  if(length >= 0.0) {
    std::size_t n =  m_procForced.size();
    for(std::size_t i=0; i<n; ++i) {
      if(procname == m_procForced[i] && r == m_regnamesForced[i] ) { 
	m_lengthForced[i] = length;
	m_weightForced[i] = wflag;
	return; 
      }
    }
    m_regnamesForced.push_back(r);
    m_procForced.push_back(procname);
    m_lengthForced.push_back(length);
    m_weightForced.push_back(wflag);
  } else {
    G4ExceptionDescription ed;
    ed << "Process: " << procname << " in region " << r
       << " : forced interacttion length= " 
       << length << " is negative - ignored"; 
    PrintWarning(ed);
  }
}

void 
G4EmExtraParameters::ActivateSecondaryBiasing(const G4String& procname,
                                              const G4String& region, 
                                              G4double factor,
                                              G4double energyLim)
{
  const G4String& r = CheckRegion(region);
  if(factor >= 0.0 && energyLim >= 0.0) {
    std::size_t n =  m_procBiasedSec.size();
    for(std::size_t i=0; i<n; ++i) {
      if(procname == m_procBiasedSec[i] && r == m_regnamesBiasedSec[i] ) { 
	m_factBiasedSec[i] = factor;
	m_elimBiasedSec[i] = energyLim;
	return; 
      }
    }
    m_regnamesBiasedSec.push_back(r);
    m_procBiasedSec.push_back(procname);
    m_factBiasedSec.push_back(factor);
    m_elimBiasedSec.push_back(energyLim);
  } else {
    G4ExceptionDescription ed;
    ed << "Process: " << procname << " in region " << r
       << " : secondary bised factor= " 
       << factor << ", Elim= " << energyLim <<  " - ignored"; 
    PrintWarning(ed);
  }
}

void G4EmExtraParameters::DefineRegParamForLoss(G4VEnergyLossProcess* ptr) const
{
  const G4RegionStore* regionStore = G4RegionStore::GetInstance();
  std::size_t n = m_regnamesSubCut.size();
  for(std::size_t i=0; i<n; ++i) { 
    const G4Region* reg = regionStore->GetRegion(m_regnamesSubCut[i], false);
    if(nullptr != reg) { ptr->ActivateSubCutoff(reg); }
  }
  n = m_procBiasedXS.size();
  for(std::size_t i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedXS[i]) {
      ptr->SetCrossSectionBiasingFactor(m_factBiasedXS[i], 
					m_weightBiasedXS[i]);
      break; 
    }
  }
  n = m_procForced.size();
  for(std::size_t i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procForced[i]) {
      ptr->ActivateForcedInteraction(m_lengthForced[i],
				     m_regnamesForced[i],
				     m_weightForced[i]);
      break; 
    }
  }
  n = m_procBiasedSec.size();
  for(std::size_t i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedSec[i]) {
      ptr->ActivateSecondaryBiasing(m_regnamesBiasedSec[i],
				    m_factBiasedSec[i], 
				    m_elimBiasedSec[i]);
      break; 
    }
  }
}

void G4EmExtraParameters::DefineRegParamForEM(G4VEmProcess* ptr) const
{
  std::size_t n = m_procBiasedXS.size();
  for(std::size_t i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedXS[i]) {
      ptr->SetCrossSectionBiasingFactor(m_factBiasedXS[i], 
					m_weightBiasedXS[i]);
      break; 
    }
  }
  n = m_procForced.size();
  for(std::size_t i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procForced[i]) {
      ptr->ActivateForcedInteraction(m_lengthForced[i],
				     m_regnamesForced[i],
				     m_weightForced[i]);
      break; 
    }
  }
  n = m_procBiasedSec.size();
  for(std::size_t i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedSec[i]) {
      ptr->ActivateSecondaryBiasing(m_regnamesBiasedSec[i],
				    m_factBiasedSec[i], 
				    m_elimBiasedSec[i]);
      break; 
    }
  }
}

G4bool G4EmExtraParameters::QuantumEntanglement()
{
  return quantumEntanglement;
}

void G4EmExtraParameters::SetQuantumEntanglement(G4bool v)
{
  quantumEntanglement = v;
}

G4bool G4EmExtraParameters::GetDirectionalSplitting() { 
  return directionalSplitting; 
}

void G4EmExtraParameters::SetDirectionalSplitting(G4bool v) 
{ 
  directionalSplitting = v; 
}

void 
G4EmExtraParameters::SetDirectionalSplittingTarget(const G4ThreeVector& v)
{ 
  directionalSplittingTarget = v; 
}

G4ThreeVector G4EmExtraParameters::GetDirectionalSplittingTarget() const
{ 
  return directionalSplittingTarget; 
}

void G4EmExtraParameters::SetDirectionalSplittingRadius(G4double r)
{ 
  directionalSplittingRadius = r; 
}

G4double G4EmExtraParameters::GetDirectionalSplittingRadius()
{ 
  return directionalSplittingRadius; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
