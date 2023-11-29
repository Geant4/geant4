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
// ClassName:  G4HadProcesses
//
// Author: 8 July 2020 V.Ivanchenko
//
// Modified:
//
// Description: access to hadronic processes via particle name or pointer
//              access or addition of hadronic cross sections
//
//----------------------------------------------------------------------------
//
#ifndef G4HadProcesses_h
#define G4HadProcesses_h 1

#include "globals.hh"
#include "G4HadronicProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionElastic.hh"

class G4VCrossSectionDataSet;

class G4HadProcesses
{
public: 

  // convert particle name to a pointer
  static const G4ParticleDefinition* FindParticle(const G4String&);

  // access to a process created before
  static G4HadronicProcess* FindInelasticProcess(const G4ParticleDefinition*);
  static G4HadronicProcess* FindInelasticProcess(const G4String&);

  static G4HadronicProcess* FindElasticProcess(const G4ParticleDefinition*);
  static G4HadronicProcess* FindElasticProcess(const G4String&);

  static G4HadronicProcess* FindCaptureProcess();

  static G4HadronicProcess* FindFissionProcess();

  // access to cross section by the component name
  static G4CrossSectionInelastic* InelasticXS(const G4String& componentName);
  static G4CrossSectionElastic* ElasticXS(const G4String& componentName);

  // add extra data set to a particle
  static G4bool AddInelasticCrossSection(const G4ParticleDefinition*, G4VCrossSectionDataSet*);
  static G4bool AddInelasticCrossSection(const G4String&, G4VCrossSectionDataSet*);

  static G4bool AddElasticCrossSection(const G4ParticleDefinition*, G4VCrossSectionDataSet*);
  static G4bool AddElasticCrossSection(const G4String&, G4VCrossSectionDataSet*);

  static G4bool AddCaptureCrossSection(G4VCrossSectionDataSet*);

  static G4bool AddFissionCrossSection(G4VCrossSectionDataSet*);

  // build neutron physics, neutron inelastic and elastic processes
  // and list of models should be defined before call to this functions
  static void BuildNeutronInelasticAndCapture(G4HadronicProcess*);
  static void BuildNeutronElastic(G4HadronicProcess*);

};

#endif

