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
// $Id: G4HadronProcessStore.cc,v 1.1 2006/10/31 11:35:02 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4HadronProcessStore
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.05.2006
//
// Modifications:
// 04.08.06 V.Ivanchenko add computation of cross sections
//
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4HadronProcessStore.hh"
#include "G4Element.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"

G4HadronProcessStore* G4HadronProcessStore::theInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronProcessStore* G4HadronProcessStore::Instance()
{
  if(0 == theInstance) {
    static G4HadronProcessStore manager;
    theInstance = &manager;
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronProcessStore::~G4HadronProcessStore()
{
  for (G4int i=0; i<n_proc; i++) {
    if( process[i] ) delete process[i];
  }
  //  for (G4int j=0; j<n_model; j++) {
  //  if( model[j] ) delete model[j];
  // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronProcessStore::G4HadronProcessStore()
{
  n_proc = 0;
  n_part = 0;
  n_model= 0;
  currentProcess  = 0;
  currentParticle = 0;
  verbose = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronProcessStore::GetInelasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      GetInelasticCrossSectionPerAtom(aParticle,kineticEnergy,elm);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronProcessStore::GetInelasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement)
{
  G4HadronicProcess* hp = FindInelasticProcess(aParticle);
  localDP.SetDefinition(const_cast<G4ParticleDefinition*>(aParticle));
  localDP.SetKineticEnergy(kineticEnergy);
  G4double cross = 0.0;
  if(hp) cross = hp->GetMicroscopicCrossSection(&localDP,
						anElement,
						STP_Temperature);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronProcessStore::GetInelasticCrossSectionPerIsotope(
    const G4ParticleDefinition *,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronProcessStore::GetElasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      GetElasticCrossSectionPerAtom(aParticle,kineticEnergy,elm);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronProcessStore::GetElasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement)
{
  G4HadronicProcess* hp = FindElasticProcess(aParticle);
  localDP.SetKineticEnergy(kineticEnergy);
  G4double cross = 0.0;
  if(hp) cross = hp->GetMicroscopicCrossSection(&localDP,
						anElement,
						STP_Temperature);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronProcessStore::GetElasticCrossSectionPerIsotope(
    const G4ParticleDefinition*,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronProcessStore::Register(G4HadronicProcess* proc, 
				    const G4ParticleDefinition* part,
				    G4HadronicInteraction* mod,
				    const G4String& name)
{
  G4int i=0;
  for(; i<n_proc; i++) {if(process[i] == proc) break;}
  G4int j=0;
  for(; j<n_part; j++) {if(particle[j] == part) break;}
  G4int k=0;
  for(; k<n_model; k++) {if(model[k] == mod) break;}
  
  if(i == n_proc || j == n_part) 
    p_map.insert(std::multimap<PD,HP>::value_type(part,proc));
 
  m_map.insert(std::multimap<HP,HI>::value_type(proc,mod));
    
  if(i == n_proc) {
    n_proc++;
    process.push_back(proc);
  }
  if(j == n_part) {
    n_part++;
    particle.push_back(part);
  }
  if(k == n_model) {
    n_model++;
    model.push_back(mod);
    modelName.push_back(name);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronProcessStore::Print(const G4ParticleDefinition* part)
{
  G4cout<<G4endl;
  G4cout << "       Hadronic Processes for " 
	 <<part->GetParticleName() << G4endl; 
  HP hp = 0;
  HI hi = 0;
  G4bool first;
  std::multimap<PD,HP,std::less<PD> >::iterator it;
  std::multimap<HP,HI,std::less<HP> >::iterator ih;
  for(it=p_map.lower_bound(part); it!=p_map.upper_bound(part); ++it) {
    if(it->first == part) {
      hp = it->second;
      G4cout << std::setw(10) << hp->GetProcessName()  
	     << "           Models:  ";
      first = true;
      for(ih=m_map.lower_bound(hp); ih!=m_map.upper_bound(hp); ++ih) {
	if(ih->first == hp) {
	  hi = ih->second;
          G4int i=0;
	  for(; i<n_model; i++) {
	    if(model[i] == hi) break;
	  }
          if(!first) G4cout << "                              ";
          first = false;
          G4cout << std::setw(10) << modelName[i] 
		 << "      Emin(GeV)= "  
		 << std::setw(5) << hi->GetMinEnergy()/GeV
		 << "      Emax(GeV)= " 
		 << std::setw(7) << hi->GetMaxEnergy()/GeV
		 << G4endl;
	}
      }
    }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronProcessStore::Dump(G4int level)
{
  if(level > 0) 
    G4cout << "=============================================================="
	   << "============================="
	   << G4endl;
  for(G4int i=0; i<n_part; i++) {
    G4String pname = (particle[i])->GetParticleName();
    G4bool yes = false;
    if(level >= 2) yes = true;
    else if(level == 1 && (pname == "proton" || 
			   pname == "neutron" ||
			   pname == "pi+" ||
			   pname == "pi-" ||
			   pname == "kaon+" ||
			   pname == "kaon-" ||
			   pname == "lambda" ||
			   pname == "anti_neutron" ||
			   pname == "anti_proton")) yes = true;
    if(yes) Print(particle[i]);
  }
  if(level > 0) 
    G4cout << "=============================================================="
	   << "============================="
	   << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronProcessStore::SetVerbose(G4int val)
{
  verbose = val;
  G4int i;
  for(i=0; i<n_proc; i++) {
    if(process[i]) process[i]->SetVerboseLevel(val);
  }
  for(i=0; i<n_model; i++) {
    if(model[i]) model[i]->SetVerboseLevel(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4int G4HadronProcessStore::GetVerbose()
{
  return verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronicProcess* G4HadronProcessStore::FindElasticProcess(
		   const G4ParticleDefinition* part)
{
  bool isNew = false;
  G4HadronicProcess* hp = 0;

  if(part != currentParticle) {
    isNew = true;
    currentParticle = part;
    localDP.SetDefinition(const_cast<G4ParticleDefinition*>(part));
  } else if(!currentProcess) {
    isNew = true;
  } else if(currentProcess->GetProcessName() == "hElastic") {
    hp = currentProcess;
  } else {
    isNew = true;
  }

  if(isNew) {
    std::multimap<PD,HP,std::less<PD> >::iterator it;
    for(it=p_map.lower_bound(part); it!=p_map.upper_bound(part); ++it) {
      if(it->first == part && (it->second)->GetProcessName() == "hElastic") {
	hp = it->second;
	break;
      }
    }  
    currentProcess = hp;
  }

  return hp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronicProcess* G4HadronProcessStore::FindInelasticProcess(
                   const G4ParticleDefinition* part)
{
  bool isNew = false;
  G4HadronicProcess* hp = 0;

  if(part != currentParticle) {
    isNew = true;
    currentParticle = part;
    localDP.SetDefinition(const_cast<G4ParticleDefinition*>(part));
  } else if(!currentProcess) {
    isNew = true;
  } else if(currentProcess->GetProcessName() == "hInelastic") {
    hp = currentProcess;
  } else {
    isNew = true;
  }

  if(isNew) {
    std::multimap<PD,HP,std::less<PD> >::iterator it;
    for(it=p_map.lower_bound(part); it!=p_map.upper_bound(part); ++it) {
      if(it->first == part && (it->second)->GetProcessName() == "hInelastic") {
	hp = it->second;
	break;
      }
    }  
    currentProcess = hp;
  }

  return hp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
