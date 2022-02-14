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

void G4PhysListUtil::InitialiseParameters()
{
  G4NistManager::Instance();
  G4EmParameters::Instance()->SetDefaults();
  G4HadronicParameters::Instance();
  G4NuclearLevelData::GetInstance();
}

G4HadronicProcess* 
G4PhysListUtil::FindInelasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = nullptr;
  if(p) {
     G4ProcessVector*  pvec = p->GetProcessManager()->GetProcessList();
     size_t n = pvec->size();
     for(size_t i=0; i<n; ++i) {
       auto proc = (*pvec)[i];
       if(proc != nullptr && fHadronInelastic == proc->GetProcessSubType()) {
	 had = static_cast<G4HadronicProcess*>(proc);
	 break;
       }
     }
  }
  return had;
}

G4HadronicProcess* 
G4PhysListUtil::FindElasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = nullptr;
  if(p) {
     G4ProcessVector*  pvec = p->GetProcessManager()->GetProcessList();
     size_t n = pvec->size();
     for(size_t i=0; i<n; ++i) {
       auto proc = (*pvec)[i];
       if(proc != nullptr && fHadronElastic == proc->GetProcessSubType()) {
	 had = static_cast<G4HadronicProcess*>(proc);
	 break;
       }
     }
  }
  return had;
}

G4HadronicProcess* 
G4PhysListUtil::FindCaptureProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = nullptr;
  if(p) {
     G4ProcessVector*  pvec = p->GetProcessManager()->GetProcessList();
     size_t n = pvec->size();
     for(size_t i=0; i<n; ++i) {
       auto proc = (*pvec)[i];
       if(proc != nullptr && fCapture == proc->GetProcessSubType()) {
	 had = static_cast<G4HadronicProcess*>(proc);
	 break;
       }
     }
  }
  return had;
}

G4HadronicProcess* 
G4PhysListUtil::FindFissionProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = nullptr;
  if(p) {
     G4ProcessVector*  pvec = p->GetProcessManager()->GetProcessList();
     size_t n = pvec->size();
     for(size_t i=0; i<n; ++i) {
       auto proc = (*pvec)[i];
       if(proc != nullptr && fFission == proc->GetProcessSubType()) {
	 had = static_cast<G4HadronicProcess*>(proc);
	 break;
       }
     }
  }
  return had;
}
