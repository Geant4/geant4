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
// $Id: G4HadronProcessStore.hh,v 1.1 2006/10/31 11:35:01 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
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

#ifndef G4HadronProcessStore_h
#define G4HadronProcessStore_h 1


#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleDefinition.hh"
#include <map>
#include <vector>

class G4Element;

class G4HadronProcessStore
{

public:

  static G4HadronProcessStore* Instance();

  ~G4HadronProcessStore();

  G4double GetInelasticCrossSectionPerVolume(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Material *material);

  G4double GetInelasticCrossSectionPerAtom(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    const G4Element *anElement);

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
    const G4Element *anElement);

  G4double GetElasticCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  void Register(G4HadronicProcess*, 
		const G4ParticleDefinition*,
		G4HadronicInteraction*,
		const G4String& name);

  void Print(const G4ParticleDefinition*);

  void Dump(G4int level);

  void SetVerbose(G4int val);

  G4int GetVerbose();

  G4HadronicProcess* FindElasticProcess(const G4ParticleDefinition*);

  G4HadronicProcess* FindInelasticProcess(const G4ParticleDefinition*);

private:

  G4HadronProcessStore();

  static G4HadronProcessStore* theInstance;

  typedef const G4ParticleDefinition* PD;
  typedef G4HadronicProcess* HP;
  typedef G4HadronicInteraction* HI;

  std::multimap<PD,HP> p_map;
  std::multimap<HP,HI> m_map;

  std::vector<G4HadronicProcess*> process;
  std::vector<G4HadronicInteraction*> model;
  std::vector<G4String> modelName;
  std::vector<PD> particle;

  HP   currentProcess;
  PD   currentParticle;

  G4DynamicParticle localDP;

  G4int n_proc;
  G4int n_model;
  G4int n_part;

  G4int verbose;

};


#endif

