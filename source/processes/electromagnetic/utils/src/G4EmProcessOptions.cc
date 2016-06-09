//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4EmProcessOptions.cc,v 1.2 2004/05/17 09:46:56 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-02 $
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
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetLossFluctuations(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSubCutoff(G4bool val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetSubCutoff(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetIntegral(G4bool val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetIntegral(val);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->SetIntegral(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMinSubRange(G4double val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetMinSubRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMinEnergy(G4double val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetMinKinEnergy(val);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->SetMinKinEnergy(val);
  }
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    G4VMultipleScattering* s = *itm;
    if(s) s->SetMinKinEnergy(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergy(G4double val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetMaxKinEnergy(val);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->SetMaxKinEnergy(val);
  }
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    G4VMultipleScattering* s = *itm;
    if(s) s->SetMaxKinEnergy(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergyForPreciseRange(G4double val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetMaxKinEnergyForPreciseRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergyForMuons(G4double val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) {
      if(abs(p->Particle()->GetPDGMass() - 105.66*MeV) < MeV)
        p->SetMaxKinEnergy(val);
    }
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) {
      if(abs(q->Particle()->GetPDGMass() - 105.66*MeV) < MeV)
        q->SetMaxKinEnergy(val);
    }
  }
  /*
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    G4VMultipleScattering* s = *itm;
    if(s) {
      if(abs(s->Particle()->GetPDGMass() - 105.66*MeV) < MeV)
        s->SetMaxKinEnergy(val);
    }
  }
  */
  theManager->SetMaxEnergyForMuons(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDEDXBinning(G4int val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetDEDXBinning(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDEDXBinningForPreciseRange(G4int val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetDEDXBinningForPreciseRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLambdaBinning(G4int val)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetLambdaBinning(val);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->SetLambdaBinning(val);
  }
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    G4VMultipleScattering* s = *itm;
    if(s) s->SetBinning(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetStepLimits(G4double v1, G4double v2)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->SetStepFunction(v1, v2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetBuildPreciseRange(G4bool val)
{
  theManager->SetBuildPreciseRange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetVerbose(G4int val, const G4String& name)
{
  G4bool all = false;
  if("all" == name) all = true;
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) {
      if(all) p->SetVerboseLevel(val);
      else if (p->GetProcessName() == name) p->SetVerboseLevel(val);
    }
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) {
      if(all) q->SetVerboseLevel(val);
      else if (q->GetProcessName() == name) q->SetVerboseLevel(val);
    }
  }
  const std::vector<G4VMultipleScattering*>& u =
        theManager->GetMultipleScatteringVector();
  std::vector<G4VMultipleScattering*>::const_iterator itm;
  for(itm = u.begin(); itm != u.end(); itm++) {
    G4VMultipleScattering* s = *itm;
    if(s) {
      if(all) s->SetVerboseLevel(val);
      else if (s->GetProcessName() == name) s->SetVerboseLevel(val);
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
    if(p) p->SetLambdaFactor(val);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->SetLambdaFactor(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::ActivateFluorescence(G4bool val, const G4Region* r)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->ActivateFluorescence(val,r);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->ActivateFluorescence(val,r);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::ActivateAugerElectronProduction(G4bool val, const G4Region* r)
{
  const std::vector<G4VEnergyLossProcess*>& v =
        theManager->GetEnergyLossProcessVector();
  std::vector<G4VEnergyLossProcess*>::const_iterator itr;
  for(itr = v.begin(); itr != v.end(); itr++) {
    G4VEnergyLossProcess* p = *itr;
    if(p) p->ActivateAugerElectronProduction(val,r);
  }
  const std::vector<G4VEmProcess*>& w =
        theManager->GetEmProcessVector();
  std::vector<G4VEmProcess*>::const_iterator itp;
  for(itp = w.begin(); itp != w.end(); itp++) {
    G4VEmProcess* q = *itp;
    if(q) q->ActivateAugerElectronProduction(val,r);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

