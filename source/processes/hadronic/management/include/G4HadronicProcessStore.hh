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
// $Id: G4HadronicProcessStore.hh 90394 2015-05-27 12:14:48Z gcosmo $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4HadronicProcessStore
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 09.05.2008
//
// Modifications:
//
//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4HadronicProcessStore_h
#define G4HadronicProcessStore_h 1


#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcessType.hh"
#include "G4ThreadLocalSingleton.hh"
#include <map>
#include <vector>
#include <iostream>

class G4Element;
class G4HadronicEPTestMessenger;

class G4HadronicProcessStore
{

friend class G4ThreadLocalSingleton<G4HadronicProcessStore>;

public:

  static G4HadronicProcessStore* Instance();

  ~G4HadronicProcessStore();

  void Clean();
  G4double GetCrossSectionPerAtom(
    const G4ParticleDefinition* particle,
    G4double kineticEnergy,
    const G4VProcess* process,
    const G4Element*  element,
    const G4Material* material=0);
      
  G4double GetCrossSectionPerVolume(
    const G4ParticleDefinition* particle,
    G4double kineticEnergy,
    const G4VProcess* process,
    const G4Material* material);
    
  G4double GetInelasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material);

  G4double GetInelasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat=0);

  G4double GetInelasticCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  G4double GetElasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material);

  G4double GetElasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat=0);

  G4double GetElasticCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  G4double GetCaptureCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material);

  G4double GetCaptureCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat=0);

  G4double GetCaptureCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  G4double GetFissionCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material);

  G4double GetFissionCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat=0);

  G4double GetFissionCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  G4double GetChargeExchangeCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material);

  G4double GetChargeExchangeCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement, const G4Material* mat=0);

  G4double GetChargeExchangeCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  // register/deregister processes following G4HadronicProcess interface
  void Register(G4HadronicProcess*); 

  void RegisterParticle(G4HadronicProcess*,
			const G4ParticleDefinition*); 

  void RegisterInteraction(G4HadronicProcess*,
			   G4HadronicInteraction*); 

  void DeRegister(G4HadronicProcess*); 

  // register/deregister processes following only G4VProcess interface
  void RegisterExtraProcess(G4VProcess*); 

  void RegisterParticleForExtraProcess(G4VProcess*,
				       const G4ParticleDefinition*); 

  void DeRegisterExtraProcess(G4VProcess*); 

  void PrintInfo(const G4ParticleDefinition*); 

  void Dump(G4int level);
  void DumpHtml();
  void PrintHtml(const G4ParticleDefinition*, std::ofstream&);
  void PrintModelHtml(const G4HadronicInteraction * model) const;

  void SetVerbose(G4int val);

  G4int GetVerbose();

  G4HadronicProcess* FindProcess(const G4ParticleDefinition*, 
				 G4HadronicProcessType subType);

  // Energy-momentum non-conservation limits and reporting
  void SetEpReportLevel(G4int level);

  void SetProcessAbsLevel(G4double absoluteLevel);

  void SetProcessRelLevel(G4double relativeLevel);

private:

  // constructor
  G4HadronicProcessStore();

  // print process info
  void Print(G4int idxProcess, G4int idxParticle);
  
  G4String HtmlFileName(const G4String &) const;

  static G4ThreadLocal G4HadronicProcessStore* instance;

  typedef const G4ParticleDefinition* PD;
  typedef G4HadronicProcess* HP;
  typedef G4HadronicInteraction* HI;

  // hadronic processes following G4HadronicProcess interface
  std::vector<G4HadronicProcess*> process;
  std::vector<G4HadronicInteraction*> model;
  std::vector<G4String> modelName;
  std::vector<PD> particle;
  std::vector<G4int> wasPrinted;

  std::multimap<PD,HP> p_map;
  std::multimap<HP,HI> m_map;

  // hadronic processes following only G4VProcess interface
  std::vector<G4VProcess*> extraProcess;
  std::multimap<PD,G4VProcess*> ep_map;

  // counters and options
  G4int n_proc;
  G4int n_model;
  G4int n_part;
  G4int n_extra;

  G4int  verbose;
  G4bool buildTableStart;

  // cache
  HP   currentProcess;
  PD   currentParticle;
  PD   theGenericIon;

  G4DynamicParticle localDP;

  G4HadronicEPTestMessenger* theEPTestMessenger;
};


#endif

