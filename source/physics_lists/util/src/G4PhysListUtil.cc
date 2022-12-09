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
//---------------------------------------------------------------------------
//
// ClassName: G4PhyslistUtil:
//     "Container" for function needed in various places  
//
// Author: 2007 Gunter Folger
//
// Modified:
// 07.10.2020 V.Ivanchenko added InitialiseParameters method
//
//----------------------------------------------------------------------------
//

#include "G4PhysListUtil.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4NistManager.hh"
#include "G4EmParameters.hh"
#include "G4HadronicParameters.hh"
#include "G4NuclearLevelData.hh"
#include "G4Neutron.hh"
#include "G4FTFTunings.hh"

void G4PhysListUtil::InitialiseParameters()
{
  G4NistManager::Instance();
  G4EmParameters::Instance();
  G4HadronicParameters::Instance();
  G4NuclearLevelData::GetInstance();
  G4FTFTunings::Instance();
}

G4VProcess* G4PhysListUtil::FindProcess(const G4ParticleDefinition* part,
                                        G4int subtype)
{
  G4VProcess* proc = nullptr;
  if(nullptr == part) { return proc; }
  G4ProcessVector* pvec = part->GetProcessManager()->GetProcessList();
  if(nullptr == pvec) { return proc; }
  G4int n = (G4int)pvec->size();
  for(G4int i=0; i<n; ++i) {
    auto ptr = (*pvec)[i];
    if(ptr != nullptr && subtype == ptr->GetProcessSubType()) {
      proc = ptr;
      break;
    }
  }
  return proc;
}

G4HadronicProcess* 
G4PhysListUtil::FindInelasticProcess(const G4ParticleDefinition* p)
{
  auto proc = FindProcess(p, fHadronInelastic);
  return dynamic_cast<G4HadronicProcess*>(proc);
}

G4HadronicProcess* 
G4PhysListUtil::FindElasticProcess(const G4ParticleDefinition* p)
{
  auto proc = FindProcess(p, fHadronElastic);
  return dynamic_cast<G4HadronicProcess*>(proc);
}

G4HadronicProcess* 
G4PhysListUtil::FindCaptureProcess(const G4ParticleDefinition* p)
{
  auto proc = FindProcess(p, fCapture);
  return dynamic_cast<G4HadronicProcess*>(proc);
}

G4HadronicProcess* 
G4PhysListUtil::FindFissionProcess(const G4ParticleDefinition* p)
{
  auto proc = FindProcess(p, fFission);
  return dynamic_cast<G4HadronicProcess*>(proc);
}

G4NeutronGeneralProcess* G4PhysListUtil::FindNeutronGeneralProcess()
{
  auto neutron = G4Neutron::Neutron();
  auto proc = dynamic_cast<G4NeutronGeneralProcess*>(FindProcess(neutron, fNeutronGeneral));
  if(nullptr == proc) {
    proc = new G4NeutronGeneralProcess();
    neutron->GetProcessManager()->AddDiscreteProcess(proc);
  }
  return proc;
}
