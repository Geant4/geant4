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
// $Id: G4HadronProcessStore.hh,v 1.1 2006-05-04 16:48:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    const G4Element *anElement);

  G4double GetInelasticCrossSectionPerIsotope(
    const G4ParticleDefinition *aParticle,
    G4double kineticEnergy,
    G4int Z, G4int A);

  G4double GetElasticCrossSectionPerVolume(
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

  void Dump();

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

  std::multimap<PD,HP,std::less<PD> > p_map;
  std::multimap<HP,HI,std::less<HP> > m_map;

  std::vector<G4HadronicProcess*> process;
  std::vector<G4HadronicInteraction*> model;
  std::vector<G4String> modelName;
  std::vector<PD> particle;

  // cash
  HP   currentProcess;
  HP   currentModel;
  PD   currentParticle;

  G4DynamicParticle localDP;

  G4int n_proc;
  G4int n_model;
  G4int n_part;

  G4int verbose;

};


#endif

