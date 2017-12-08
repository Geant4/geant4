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
// $Id: G4HadronicProcessStore.cc 105940 2017-09-01 07:33:44Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4HadronicProcessStore
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 09.05.2008
//
// Modifications:
// 23.01.2009 V.Ivanchenko add destruction of processes
//
// Class Description:
// Singleton to store hadronic processes, to provide access to processes
// and to printout information about processes
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4HadronicProcessStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Element.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4ParticleTable.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4HadronicEPTestMessenger.hh"
#ifdef G4USE_STD11
  #include <algorithm>   //transform
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4ThreadLocal G4HadronicProcessStore* G4HadronicProcessStore::instance = 0;

G4HadronicProcessStore* G4HadronicProcessStore::Instance()
{
  if(!instance) {
    static G4ThreadLocalSingleton<G4HadronicProcessStore> inst;
    instance = inst.Instance();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronicProcessStore::~G4HadronicProcessStore()
{
  Clean();
  delete theEPTestMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::Clean()
{
  G4int i;
  //std::cout << "G4HadronicProcessStore::Clean() Nproc= " << n_proc
  //	    << "  Nextra= " << n_extra << std::endl;
  for (i=0; i<n_proc; ++i) {
    if( process[i] ) {
      //G4cout << "G4HadronicProcessStore::Clean() delete hadronic "
      //  << i << "  " <<  process[i]->GetProcessName() << G4endl;
      delete process[i];
    }
  }
  for(i=0; i<n_extra; ++i) {
    if(extraProcess[i]) {
        // G4cout << "G4HadronicProcessStore::Clean() delete extra proc "
        //<< i << "  " << extraProcess[i]->GetProcessName() << G4endl;
        delete extraProcess[i];
        extraProcess[i] = 0;
    }
  }
  //std::cout << "G4HadronicProcessStore::Clean() done" << std::endl;
  n_extra = 0;
  n_proc = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronicProcessStore::G4HadronicProcessStore()
{
  n_proc = 0;
  n_part = 0;
  n_model= 0;
  n_extra= 0;
  currentProcess  = 0;
  currentParticle = 0;
  theGenericIon = 
    G4ParticleTable::GetParticleTable()->FindParticle("GenericIon");
  verbose = 1;
  buildTableStart = true;
  theEPTestMessenger = new G4HadronicEPTestMessenger(this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetCrossSectionPerAtom(
                                 const G4ParticleDefinition* part,
                                 G4double energy,
                                 const G4VProcess* proc,
                                 const G4Element*  element,
				 const G4Material* material)
{
  G4double cross = 0.;    
  G4int subType = proc->GetProcessSubType();      
  if (subType == fHadronElastic)   
    cross = GetElasticCrossSectionPerAtom(part,energy,element,material);
  else if (subType == fHadronInelastic)   
    cross = GetInelasticCrossSectionPerAtom(part,energy,element,material);
  else if (subType == fCapture)   
    cross = GetCaptureCrossSectionPerAtom(part,energy,element,material);      
  else if (subType == fFission)   
    cross = GetFissionCrossSectionPerAtom(part,energy,element,material); 
  else if (subType == fChargeExchange)   
    cross = GetChargeExchangeCrossSectionPerAtom(part,energy,element,material);
  return cross;
}  	   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetCrossSectionPerVolume(
                                 const G4ParticleDefinition* part,
                                 G4double energy,
                                 const G4VProcess* proc,
                                 const G4Material* material)
{
  G4double cross = 0.;    
  G4int subType = proc->GetProcessSubType();      
  if (subType == fHadronElastic)   
    cross = GetElasticCrossSectionPerVolume(part,energy,material);
  else if (subType == fHadronInelastic)   
    cross = GetInelasticCrossSectionPerVolume(part,energy,material);
  else if (subType == fCapture)   
    cross = GetCaptureCrossSectionPerVolume(part,energy,material);      
  else if (subType == fFission)   
    cross = GetFissionCrossSectionPerVolume(part,energy,material); 
  else if (subType == fChargeExchange)   
    cross = GetChargeExchangeCrossSectionPerVolume(part,energy,material);
  return cross;
}  	   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetElasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; ++i) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      GetElasticCrossSectionPerAtom(aParticle,kineticEnergy,elm,material);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetElasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat)
{
  G4HadronicProcess* hp = FindProcess(aParticle, fHadronElastic);
  G4double cross = 0.0;
  localDP.SetKineticEnergy(kineticEnergy);
  if(hp) {
    cross = hp->GetElementCrossSection(&localDP,anElement,mat);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetElasticCrossSectionPerIsotope(
    const G4ParticleDefinition*,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetInelasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; ++i) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      GetInelasticCrossSectionPerAtom(aParticle,kineticEnergy,elm,material);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetInelasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat)
{
  G4HadronicProcess* hp = FindProcess(aParticle, fHadronInelastic);
  localDP.SetKineticEnergy(kineticEnergy);
  G4double cross = 0.0;
  if(hp) { 
    cross = hp->GetElementCrossSection(&localDP,anElement,mat);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetInelasticCrossSectionPerIsotope(
    const G4ParticleDefinition *,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetCaptureCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; ++i) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      GetCaptureCrossSectionPerAtom(aParticle,kineticEnergy,elm,material);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetCaptureCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat)
{
  G4HadronicProcess* hp = FindProcess(aParticle, fCapture);
  localDP.SetKineticEnergy(kineticEnergy);
  G4double cross = 0.0;
  if(hp) {
    cross = hp->GetElementCrossSection(&localDP,anElement,mat);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetCaptureCrossSectionPerIsotope(
    const G4ParticleDefinition *,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetFissionCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
      GetFissionCrossSectionPerAtom(aParticle,kineticEnergy,elm,material);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetFissionCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat)
{
  G4HadronicProcess* hp = FindProcess(aParticle, fFission);
  localDP.SetKineticEnergy(kineticEnergy);
  G4double cross = 0.0;
  if(hp) {
    cross = hp->GetElementCrossSection(&localDP,anElement,mat);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetFissionCrossSectionPerIsotope(
    const G4ParticleDefinition *,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetChargeExchangeCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; ++i) {
    const G4Element* elm = (*theElementVector)[i];
    cross += theAtomNumDensityVector[i]*
    GetChargeExchangeCrossSectionPerAtom(aParticle,kineticEnergy,elm,material);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetChargeExchangeCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat)
{
  G4HadronicProcess* hp = FindProcess(aParticle, fChargeExchange);
  localDP.SetKineticEnergy(kineticEnergy);
  G4double cross = 0.0;
  if(hp) {
    cross = hp->GetElementCrossSection(&localDP,anElement,mat);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4HadronicProcessStore::GetChargeExchangeCrossSectionPerIsotope(
    const G4ParticleDefinition *,
    G4double,
    G4int, G4int)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::Register(G4HadronicProcess* proc) 
{ 
  for(G4int i=0; i<n_proc; ++i) {
    if(process[i] == proc) { return; }
  }
  if(1 < verbose) {
    G4cout << "G4HadronicProcessStore::Register hadronic " << n_proc
	   << "  " << proc->GetProcessName() << G4endl;
  }
  ++n_proc;
  process.push_back(proc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::RegisterParticle(G4HadronicProcess* proc, 
					      const G4ParticleDefinition* part) 
{ 
  G4int i=0;
  for(; i<n_proc; ++i) {if(process[i] == proc) break;}
  G4int j=0;
  for(; j<n_part; ++j) {if(particle[j] == part) break;}

  if(1 < verbose) {
    G4cout << "G4HadronicProcessStore::RegisterParticle " 
	   << part->GetParticleName()
	   << " for  " << proc->GetProcessName() << G4endl;
  }
  if(j == n_part) {
    ++n_part;
    particle.push_back(part);
    wasPrinted.push_back(0);
  }
  
  // the pair should be added?
  if(i < n_proc) {
    std::multimap<PD,HP,std::less<PD> >::iterator it;
    for(it=p_map.lower_bound(part); it!=p_map.upper_bound(part); ++it) {
      if(it->first == part) {
	HP process2 = (it->second);
	if(proc == process2) { return; }
      }
    }
  }
  
  p_map.insert(std::multimap<PD,HP>::value_type(part,proc));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::RegisterInteraction(G4HadronicProcess* proc,
						 G4HadronicInteraction* mod)
{
  G4int i=0;
  for(; i<n_proc; ++i) {if(process[i] == proc) { break; }}
  G4int k=0;
  for(; k<n_model; ++k) {if(model[k] == mod) { break; }}
   
  m_map.insert(std::multimap<HP,HI>::value_type(proc,mod));
    
  if(k == n_model) {
    ++n_model;
    model.push_back(mod);
    modelName.push_back(mod->GetModelName());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::DeRegister(G4HadronicProcess* proc)
{
  for(G4int i=0; i<n_proc; ++i) {
    if(process[i] == proc) {
      process[i] = 0;
      DeRegisterExtraProcess((G4VProcess*)proc);      
      return;
    }
  }
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::RegisterExtraProcess(G4VProcess* proc)
{
  for(G4int i=0; i<n_extra; ++i) {
    if(extraProcess[i] == proc) { return; }
  }
  G4HadronicProcess* hproc = reinterpret_cast<G4HadronicProcess*>(proc);
  if(hproc) {
    for(G4int i=0; i<n_proc; ++i) {
      if(process[i] == hproc) { return; }
    }
  }
  if(1 < verbose) {
    G4cout << "Extra Process: " << n_extra 
	   << "  " <<  proc->GetProcessName() << G4endl;
  }
  ++n_extra;
  extraProcess.push_back(proc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::RegisterParticleForExtraProcess(
                             G4VProcess* proc,
			     const G4ParticleDefinition* part)
{
  G4int i=0;
  for(; i<n_extra; ++i) { if(extraProcess[i] == proc) { break; } }
  G4int j=0;
  for(; j<n_part; ++j) { if(particle[j] == part) { break; } }

  if(j == n_part) {
    ++n_part;
    particle.push_back(part);
    wasPrinted.push_back(0);
  }
  
  // the pair should be added?
  if(i < n_extra) {
    std::multimap<PD,G4VProcess*,std::less<PD> >::iterator it;
    for(it=ep_map.lower_bound(part); it!=ep_map.upper_bound(part); ++it) {
      if(it->first == part) {
	G4VProcess* process2 = (it->second);
	if(proc == process2) { return; }
      }
    }
  }
  
  ep_map.insert(std::multimap<PD,G4VProcess*>::value_type(part,proc));
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::DeRegisterExtraProcess(G4VProcess* proc)
{
  for(G4int i=0; i<n_extra; ++i) {
    if(extraProcess[i] == proc) {
      extraProcess[i] = 0;
      if(1 < verbose) {
	G4cout << "Extra Process: " << i << "  " 
	       <<proc->GetProcessName()<< " is deregisted " << G4endl;
      }
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::PrintInfo(const G4ParticleDefinition* part) 
{
  // Trigger particle/process/model printout only when last particle is 
  // registered
  if(buildTableStart && part == particle[n_part - 1]) {
    buildTableStart = false;
    Dump(verbose);
    if (getenv("G4PhysListDocDir") ) DumpHtml();
    G4HadronicInteractionRegistry::Instance()->InitialiseModels();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::DumpHtml()
{
  // Automatic generation of html documentation page for physics lists
  // List processes, models and cross sections for the most important
  // particles in descending order of importance

  char* dirName = getenv("G4PhysListDocDir");
  char* physListName = getenv("G4PhysListName");
  if (dirName && physListName) {

    // Open output file with path name
    G4String pathName = G4String(dirName) + "/" + G4String(physListName) + ".html";
    std::ofstream outFile;
    outFile.open(pathName);

    // Write physics list summary file
    outFile << "<html>\n";
    outFile << "<head>\n";
    outFile << "<title>Physics List Summary</title>\n";
    outFile << "</head>\n";
    outFile << "<body>\n";
    outFile << "<h2> Summary of Hadronic Processes, Models and Cross Sections for Physics List "
            << G4String(physListName) << "</h2>\n";
    outFile << "<ul>\n";

    PrintHtml(G4Proton::Proton(), outFile);
    PrintHtml(G4Neutron::Neutron(), outFile);
    PrintHtml(G4PionPlus::PionPlus(), outFile); 
    PrintHtml(G4PionMinus::PionMinus(), outFile);
    PrintHtml(G4Gamma::Gamma(), outFile);
    PrintHtml(G4Electron::Electron(), outFile);
//    PrintHtml(G4MuonMinus::MuonMinus(), outFile);
    PrintHtml(G4Positron::Positron(), outFile);
    PrintHtml(G4KaonPlus::KaonPlus(), outFile);
    PrintHtml(G4KaonMinus::KaonMinus(), outFile);
    PrintHtml(G4Lambda::Lambda(), outFile);
    PrintHtml(G4Alpha::Alpha(), outFile);
    PrintHtml(G4GenericIon::GenericIon(), outFile);

    outFile << "</ul>\n";
    outFile << "</body>\n";
    outFile << "</html>\n";
    outFile.close();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::PrintHtml(const G4ParticleDefinition* theParticle,
                                       std::ofstream& outFile)
{
  // Automatic generation of html documentation page for physics lists
  // List processes for the most important particles in descending order
  // of importance
 
  outFile << "<br> <li><h2><font color=\" ff0000 \">" 
          << theParticle->GetParticleName() << "</font></h2></li>\n";

  typedef std::multimap<PD,HP,std::less<PD> > PDHPmap;
  typedef std::multimap<HP,HI,std::less<HP> > HPHImap;

  std::pair<PDHPmap::iterator, PDHPmap::iterator> itpart =
                        p_map.equal_range(theParticle);

  // Loop over processes assigned to particle

  G4HadronicProcess* theProcess;
  for (PDHPmap::iterator it = itpart.first; it != itpart.second; ++it) {
    theProcess = (*it).second;
    //  description is inline
    //outFile << "<br> &nbsp;&nbsp; <b><font color=\" 0000ff \">process : <a href=\""
    //        << theProcess->GetProcessName() << ".html\"> "
    //        << theProcess->GetProcessName() << "</a></font></b>\n";
    outFile << "<br> &nbsp;&nbsp; <b><font color=\" 0000ff \">process : "
            << theProcess->GetProcessName() << "</font></b>\n";
    outFile << "<ul>\n";
    outFile << "  <li>";
   theProcess->ProcessDescription(outFile);
   outFile << "  <li><b><font color=\" 00AA00 \">models : </font></b>\n";
    // Loop over models assigned to process
    std::pair<HPHImap::iterator, HPHImap::iterator> itmod =
                        m_map.equal_range(theProcess);

    outFile << "    <ul>\n";
	 G4String physListName(getenv("G4PhysListName"));

    for (HPHImap::iterator jt = itmod.first; jt != itmod.second; ++jt) {
      outFile << "    <li><b><a href=\"" << physListName << "_" 
		        << HtmlFileName((*jt).second->GetModelName()) << "\"> "
              << (*jt).second->GetModelName() << "</a>" 
              << " from " << (*jt).second->GetMinEnergy()/GeV
              << " GeV to " << (*jt).second->GetMaxEnergy()/GeV
              << " GeV </b></li>\n";

      // Print ModelDescription, ignore that we overwrite files n-times.
      PrintModelHtml((*jt).second);

    }
    outFile << "    </ul>\n";
    outFile << "  </li>\n";

    // List cross sections assigned to process
    outFile << "  <li><b><font color=\" 00AA00 \">cross sections : </font></b>\n";
    outFile << "    <ul>\n";
    theProcess->GetCrossSectionDataStore()->DumpHtml(*theParticle, outFile);
    //        << " \n";
    outFile << "    </ul>\n";

    outFile << "  </li>\n";
    outFile << "</ul>\n";

  }

  // Loop over extra (G4VProcess) processes

  std::multimap<PD,G4VProcess*,std::less<PD> >::iterator itp;
  for (itp=ep_map.lower_bound(theParticle); itp!=ep_map.upper_bound(theParticle); ++itp) {
    if (itp->first == theParticle) {
      G4VProcess* proc = (itp->second);
      outFile << "<br> &nbsp;&nbsp; <b><font color=\" 0000ff \">process : "
              << proc->GetProcessName() << "</font></b>\n";
      outFile << "<ul>\n";
      outFile << "  <li>";
      proc->ProcessDescription(outFile);
      outFile << "  </li>\n";
      outFile << "</ul>\n";
    }
  }

} // PrintHtml for particle

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4HadronicProcessStore::PrintModelHtml(const G4HadronicInteraction * mod) const
{
	G4String dirName(getenv("G4PhysListDocDir"));
	G4String physListName(getenv("G4PhysListName"));
	G4String pathName = dirName + "/" + physListName + "_" + HtmlFileName(mod->GetModelName());
	std::ofstream outModel;
	outModel.open(pathName);
	outModel << "<html>\n";
	outModel << "<head>\n";
	outModel << "<title>Description of " << mod->GetModelName() 
		 << "</title>\n";
	outModel << "</head>\n";
	outModel << "<body>\n";

	mod->ModelDescription(outModel);

	outModel << "</body>\n";
	outModel << "</html>\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//private 
G4String G4HadronicProcessStore::HtmlFileName(const G4String & in) const
{
   G4String str(in);
    // replace blanks by _  C++11 version:
#ifdef G4USE_STD11
	std::transform(str.begin(), str.end(), str.begin(), [](char ch) {
     return ch == ' ' ? '_' : ch;
   });
#else	
	  // and now in ancient language
	   for(std::string::iterator it = str.begin(); it != str.end(); ++it) {
        if(*it == ' ') *it = '_';
      }
#endif
   str=str + ".html";		
   return str;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::Dump(G4int level)
{
 if (level == 0) return;
  
 G4cout 
   << "\n====================================================================\n"
   << std::setw(60) << "HADRONIC PROCESSES SUMMARY (verbose level " << level
   << ")" << G4endl;
  
 for (G4int i=0; i<n_part; ++i) {
    PD part = particle[i];
    G4String pname = part->GetParticleName();
    G4bool yes = false;

    if (level == 1 && (pname == "proton" || 
		       pname == "neutron" ||
                       pname == "deuteron" ||
                       pname == "triton" ||
                       pname == "He3" ||
                       pname == "alpha" ||
		       pname == "pi+" ||
		       pname == "pi-" ||
                       pname == "gamma" ||
                       pname == "e+" ||
                       pname == "e-" ||
                       pname == "mu+" ||
                       pname == "mu-" ||
		       pname == "kaon+" ||
		       pname == "kaon-" ||
		       pname == "lambda" ||
		       pname == "GenericIon" ||
		       pname == "anti_neutron" ||
		       pname == "anti_proton" ||
                       pname == "anti_deuteron" ||
                       pname == "anti_triton" ||
                       pname == "anti_He3" ||
                       pname == "anti_alpha")) yes = true;
    if (level > 1) yes = true;			   
    if (yes) {
      // main processes
      std::multimap<PD,HP,std::less<PD> >::iterator it;

      for (it=p_map.lower_bound(part); it!=p_map.upper_bound(part); ++it) {
	if (it->first == part) {
	  HP proc = (it->second);
	  G4int j=0;
	  for (; j<n_proc; ++j) {
	    if (process[j] == proc) { Print(j, i); }
	  }
	}
      }
      
      // extra processes
      std::multimap<PD,G4VProcess*,std::less<PD> >::iterator itp;
      for(itp=ep_map.lower_bound(part); itp!=ep_map.upper_bound(part); ++itp) {
	if(itp->first == part) {
	  G4VProcess* proc = (itp->second);
	  if (wasPrinted[i] == 0) {
            G4cout << "\n---------------------------------------------------\n"
           << std::setw(50) << "Hadronic Processes for " 
	   << part->GetParticleName() << "\n";	  
	    wasPrinted[i] = 1;
	  }
	  G4cout << "\n  Process: " << proc->GetProcessName() << G4endl;
	}
      }
    }
  }

  G4cout << "\n================================================================"
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::Print(G4int idxProc, G4int idxPart)
{
  G4HadronicProcess* proc = process[idxProc];
  const G4ParticleDefinition* part = particle[idxPart];
  if (wasPrinted[idxPart] == 0) {
    G4cout << "\n---------------------------------------------------\n"  
           << std::setw(50) << "Hadronic Processes for " 
	   << part->GetParticleName() << "\n";
    wasPrinted[idxPart] = 1;	   
  }
  
  G4cout << "\n  Process: " << proc->GetProcessName();

  // Append the string "/n" (i.e. "per nucleon") on the kinetic energy of ions.
  G4String stringEnergyPerNucleon = "";
  if ( part &&
       ( part == G4GenericIon::Definition() || 
         std::abs( part->GetBaryonNumber() ) > 1 ) ) {
    stringEnergyPerNucleon = "/n";
  }

  HI hi = 0;
  std::multimap<HP,HI,std::less<HP> >::iterator ih;
  for(ih=m_map.lower_bound(proc); ih!=m_map.upper_bound(proc); ++ih) {
    if(ih->first == proc) {
      hi = ih->second;
      G4int i=0;
      for(; i<n_model; ++i) {
	if(model[i] == hi) { break; }
      }
      G4cout << "\n        Model: " << std::setw(25) << modelName[i] << ": "  
	     << G4BestUnit(hi->GetMinEnergy(), "Energy") << stringEnergyPerNucleon
	     << " ---> " 
	     << G4BestUnit(hi->GetMaxEnergy(), "Energy") << stringEnergyPerNucleon;
    }
  }
  G4cout << G4endl;
  
  G4CrossSectionDataStore* csds = proc->GetCrossSectionDataStore();
  csds->DumpPhysicsTable(*part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::SetVerbose(G4int val)
{
  verbose = val;
  G4int i;
  for(i=0; i<n_proc; ++i) {
    if(process[i]) { process[i]->SetVerboseLevel(val); }
  }
  for(i=0; i<n_model; ++i) {
    if(model[i]) { model[i]->SetVerboseLevel(val); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4int G4HadronicProcessStore::GetVerbose()
{
  return verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4HadronicProcess* G4HadronicProcessStore::FindProcess(
   const G4ParticleDefinition* part, G4HadronicProcessType subType)
{
  bool isNew = false;
  G4HadronicProcess* hp = 0;
  localDP.SetDefinition(part);

  if(part != currentParticle) {
    const G4ParticleDefinition* p = part;
    if(p->GetBaryonNumber() > 4 && p->GetParticleType() == "nucleus") {
      p = theGenericIon;
    }
    if(p !=  currentParticle) { 
      isNew = true;
      currentParticle = p;
    }
  }
  if(!isNew) { 
    if(!currentProcess) {
      isNew = true;
    } else if(subType == currentProcess->GetProcessSubType()) {
      hp = currentProcess;
    } else {
      isNew = true;
    }
  }
  if(isNew) {
    std::multimap<PD,HP,std::less<PD> >::iterator it;
    for(it=p_map.lower_bound(currentParticle); 
	it!=p_map.upper_bound(currentParticle); ++it) {
      if(it->first == currentParticle && 
	 subType == (it->second)->GetProcessSubType()) {
	hp = it->second;
	break;
      }
    }  
    currentProcess = hp;
  }
  return hp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::SetEpReportLevel(G4int level)
{
  G4cout << " Setting energy/momentum report level to " << level 
         << " for " << process.size() << " hadronic processes " << G4endl;
  for (G4int i = 0; i < G4int(process.size()); ++i) {
    process[i]->SetEpReportLevel(level);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::SetProcessAbsLevel(G4double abslevel)
{
  G4cout << " Setting absolute energy/momentum test level to " << abslevel 
	 << G4endl;
  G4double rellevel = 0.0;
  G4HadronicProcess* theProcess = 0;
  for (G4int i = 0; i < G4int(process.size()); ++i) {
    theProcess = process[i];
    rellevel = theProcess->GetEnergyMomentumCheckLevels().first;
    theProcess->SetEnergyMomentumCheckLevels(rellevel, abslevel);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4HadronicProcessStore::SetProcessRelLevel(G4double rellevel)
{
  G4cout << " Setting relative energy/momentum test level to " << rellevel 
	 << G4endl;
  G4double abslevel = 0.0;
  G4HadronicProcess* theProcess = 0;
  for (G4int i = 0; i < G4int(process.size()); ++i) {
    theProcess = process[i];
    abslevel = theProcess->GetEnergyMomentumCheckLevels().second;
    theProcess->SetEnergyMomentumCheckLevels(rellevel, abslevel);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
