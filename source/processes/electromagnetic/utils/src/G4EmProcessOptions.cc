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
// $Id: G4EmProcessOptions.cc,v 1.30 2010-11-23 19:01:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmProcessOptions
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 27.02.2004
//
// Modifications:
// 30-06-04 G4EmProcess is pure discrete (V.Ivanchenko)
// 24-03-05 Add ApplyCuts and RandomStep (V.Ivanchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 10-05-06 Add command MscStepLimit to G4LossTableManager (V.Ivantchenko) 
// 22-05-06 Add SetBremsstrahlungTh (V.Ivanchenko)
// 12-02-07 Add SetSkin, SetLinearLossLimit (V.Ivanchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmProcessOptions.hh"
#include "G4LossTableManager.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VMultipleScattering.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4VAtomDeexcitation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmProcessOptions::G4EmProcessOptions()
{
  theManager = G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmProcessOptions::~G4EmProcessOptions()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLossFluctuations(G4bool val)
{
  theManager->SetLossFluctuations(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSubCutoff(G4bool val, const G4Region* r)
{
  theManager->SetSubCutoff(val, r);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetIntegral(G4bool val)
{
  theManager->SetIntegral(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMinSubRange(G4double val)
{
  theManager->SetMinSubRange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMinEnergy(G4double val)
{
  theManager->SetMinEnergy(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergy(G4double val)
{
  theManager->SetMaxEnergy(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergyForCSDARange(G4double val)
{
  theManager->SetMaxEnergyForCSDARange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergyForMuons(G4double val)
{
  theManager->SetMaxEnergyForMuons(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDEDXBinning(G4int val)
{
  theManager->SetDEDXBinning(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDEDXBinningForCSDARange(G4int val)
{
  theManager->SetDEDXBinningForCSDARange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLambdaBinning(G4int val)
{
  theManager->SetLambdaBinning(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetStepFunction(G4double v1, G4double v2)
{
  theManager->SetStepFunction(v1, v2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetRandomStep(G4bool val)
{
  theManager->SetRandomStep(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetApplyCuts(G4bool val)
{
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) { q->SetApplyCuts(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetBuildCSDARange(G4bool val)
{
  theManager->SetBuildCSDARange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetVerbose(G4int val, const G4String& name)
{
  G4bool all = false;
  if("all" == name) { all = true; }
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();

  if(all) { 
    theManager->SetVerbose(val);
    return;
  }

  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); ++itr) {
    G4VEnergyLossProcess* p = *itr;
    if(p) {
      if (p->GetProcessName() == name) { p->SetVerboseLevel(val); }
    }
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) {
      if (q->GetProcessName() == name) { q->SetVerboseLevel(val); }
    }
  }
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    G4VMultipleScattering* s = *itm;
    if(s) {
      if (s->GetProcessName() == name) { s->SetVerboseLevel(val); }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLambdaFactor(G4double val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) { p->SetLambdaFactor(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::ActivateDeexcitation(const G4String& pname, 
					      G4bool val, 
					      const G4String& reg)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  const G4Region* r = 0;
  if(reg == "" || reg == "World") {
    r = regionStore->GetRegion("DefaultRegionForTheWorld", false);
  } else {
    r = regionStore->GetRegion(reg, false);
  }
  if(!r) {
    G4cout << "G4EmProcessOptions::ActivateDeexcitation ERROR: G4Region <"
	   << reg << "> not found, the command ignored" << G4endl;
    return;
  }

  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) {
      if(pname == p->GetProcessName()) { p->ActivateDeexcitation(val,r); }
    }
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) {
      if(pname == q->GetProcessName()) { q->ActivateDeexcitation(val,r); }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDeexcitationActive(G4bool val)
{
  G4VAtomDeexcitation* ad = theManager-> AtomDeexcitation();
  if(ad) { ad->SetActive(val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::SetDeexcitationActiveRegion(const G4String& rname, 
						G4bool valDeexcitation,
						G4bool valAuger,
						G4bool valPIXE)
{
  G4VAtomDeexcitation* ad = theManager-> AtomDeexcitation();
  if(ad) { 
    ad->SetDeexcitationActiveRegion(rname, valDeexcitation,
				    valAuger,valPIXE); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetAugerActive(G4bool val)
{
  G4VAtomDeexcitation* ad = theManager-> AtomDeexcitation();
  if(ad) { ad->SetAugerActive(val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetPIXEActive(G4bool val)
{
  G4VAtomDeexcitation* ad = theManager-> AtomDeexcitation();
  if(ad) { ad->SetPIXEActive(val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetPIXECrossSectionModel(const G4String& mname)
{
  G4VAtomDeexcitation* ad = theManager-> AtomDeexcitation();
  if(ad) { ad->SetPIXECrossSectionModel(mname); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscStepLimitation(G4MscStepLimitType val)
{
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    if(*itm) (*itm)->SetStepLimitType(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscLateralDisplacement(G4bool val)
{
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    if(*itm) { (*itm)->SetLateralDisplasmentFlag(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSkin(G4double val)
{
  if(val < 0.0) return;
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    if(*itm) { (*itm)->SetSkin(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscRangeFactor(G4double val)
{
  if(val < 0.0) return;
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    if(*itm) { (*itm)->SetRangeFactor(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscGeomFactor(G4double val)
{
  if(val < 0.0) { return; }
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    if(*itm) { (*itm)->SetGeomFactor(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetPolarAngleLimit(G4double val)
{
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    if(*itm) { (*itm)->SetPolarAngleLimit(val); }
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) { q->SetPolarAngleLimit(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetFactorForAngleLimit(G4double val)
{
  theManager->SetFactorForAngleLimit(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLPMFlag(G4bool val)
{
  theManager->SetLPMFlag(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSplineFlag(G4bool val)
{
  theManager->SetSplineFlag(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLinearLossLimit(G4double val)
{
  theManager->SetLinearLossLimit(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetBremsstrahlungTh(G4double val)
{
  theManager->SetBremsstrahlungTh(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

