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
//---------------------------------------------------------------------------
//
// ClassName:  G4HadProcesses
//
// Author: 8 July 2020 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadProcesses.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysListUtil.hh"
#include "G4ParticleTable.hh"
#include "G4Neutron.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysicsListHelper.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4NeutronCaptureXS.hh"

const G4ParticleDefinition* G4HadProcesses::FindParticle(const G4String& pname)
{
  return G4ParticleTable::GetParticleTable()->FindParticle(pname);
}

G4HadronicProcess* G4HadProcesses::FindInelasticProcess(const G4ParticleDefinition* ptr)
{
  return G4PhysListUtil::FindInelasticProcess(ptr);
}

G4HadronicProcess* G4HadProcesses::FindInelasticProcess(const G4String& pname)
{
  return FindInelasticProcess( FindParticle(pname) );
}

G4HadronicProcess* G4HadProcesses::FindElasticProcess(const G4ParticleDefinition* ptr)
{
  return G4PhysListUtil::FindElasticProcess(ptr);
}

G4HadronicProcess* G4HadProcesses::FindElasticProcess(const G4String& pname)
{
  return FindElasticProcess( FindParticle(pname) );
}

G4HadronicProcess* G4HadProcesses::FindCaptureProcess()
{
  return G4PhysListUtil::FindCaptureProcess(G4Neutron::Neutron());
}

G4HadronicProcess* G4HadProcesses::FindFissionProcess()
{
  return G4PhysListUtil::FindFissionProcess(G4Neutron::Neutron());
}

G4CrossSectionInelastic* G4HadProcesses::InelasticXS(const G4String& compName)
{
  G4CrossSectionInelastic* xs = nullptr;
  auto comp = G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection(compName);
  if( comp != nullptr ) {
    xs = new G4CrossSectionInelastic(comp);
  } else if( "Glauber-Gribov" == compName ) {
    xs = new G4CrossSectionInelastic(new G4ComponentGGHadronNucleusXsc());
  } else if( "Glauber-Gribov Nucl-nucl" == compName ) {
    xs = new G4CrossSectionInelastic(new G4ComponentGGNuclNuclXsc());
  } else if( "AntiAGlauber" == compName ) {
    xs = new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS());
  }
  return xs;
}

G4CrossSectionElastic* G4HadProcesses::ElasticXS(const G4String& compName)
{
  G4CrossSectionElastic* xs = nullptr;
  auto comp = G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection(compName);
  if( comp != nullptr ) {
    xs = new G4CrossSectionElastic(comp);
  } else if( "Glauber-Gribov" == compName ) {
    xs = new G4CrossSectionElastic(new G4ComponentGGHadronNucleusXsc());
  } else if( "Glauber-Gribov Nucl-nucl" == compName ) {
    xs = new G4CrossSectionElastic(new G4ComponentGGNuclNuclXsc());
  } else if( "AntiAGlauber" == compName ) {
    xs = new G4CrossSectionElastic(new G4ComponentAntiNuclNuclearXS());
  }
  return xs;
}

G4bool G4HadProcesses::AddInelasticCrossSection(const G4ParticleDefinition* ptr,
                                                G4VCrossSectionDataSet* xs)
{
  G4bool isOK(false);
  if( ptr != nullptr ) {
    G4HadronicProcess* had = FindInelasticProcess( ptr );
    if( had != nullptr ) {
      isOK = true;
      had->AddDataSet( xs );
    }
  }
  return isOK;
}

G4bool G4HadProcesses::AddInelasticCrossSection(const G4String& pname, G4VCrossSectionDataSet* xs)
{
  return AddInelasticCrossSection( FindParticle(pname), xs );
}

G4bool G4HadProcesses::AddElasticCrossSection(const G4ParticleDefinition* ptr, G4VCrossSectionDataSet* xs)
{
  G4bool isOK(false);
  if( ptr != nullptr ) {
    G4HadronicProcess* had = FindElasticProcess( ptr );
    if( had != nullptr ) {
      isOK = true;
      had->AddDataSet( xs );
    }
  }
  return isOK;
}

G4bool G4HadProcesses::AddElasticCrossSection(const G4String& pname, G4VCrossSectionDataSet* xs)
{
  return AddElasticCrossSection( FindParticle(pname), xs );
}

G4bool G4HadProcesses::AddCaptureCrossSection(G4VCrossSectionDataSet* xs)
{
  G4bool isOK(false);
  G4HadronicProcess* had = FindCaptureProcess();
  if( had != nullptr ) {
    isOK = true;
    had->AddDataSet( xs );
  }
  return isOK;
}

G4bool G4HadProcesses::AddFissionCrossSection(G4VCrossSectionDataSet* xs)
{
  G4bool isOK(false);
  G4HadronicProcess* had = FindFissionProcess();
  if( had != nullptr ) {
    isOK = true;
    had->AddDataSet( xs );
  }
  return isOK;
}

void G4HadProcesses::BuildNeutronInelasticAndCapture(G4HadronicProcess* nInel)
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useNeutronGeneral = param->EnableNeutronGeneralProcess();

  G4HadronicProcess* nCap = new G4NeutronCaptureProcess("nCapture");
  nCap->RegisterMe(new G4NeutronRadCapture());

  if(useNeutronGeneral) {
    auto nGen = G4PhysListUtil::FindNeutronGeneralProcess();
    nGen->SetInelasticProcess(nInel);
    nGen->SetCaptureProcess(nCap);
  } else {
    auto neutron = G4Neutron::Neutron();
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    nInel->AddDataSet(new G4NeutronInelasticXS());
    ph->RegisterProcess(nInel, neutron);
    ph->RegisterProcess(nCap, neutron);
    if( param->ApplyFactorXS() ) {
      nInel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
    }
  }
}

void G4HadProcesses::BuildNeutronElastic(G4HadronicProcess* nEl)
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useNeutronGeneral = param->EnableNeutronGeneralProcess();

  if(useNeutronGeneral) {
    auto nGen = G4PhysListUtil::FindNeutronGeneralProcess();
    nGen->SetElasticProcess(nEl);
  } else {
    auto neutron = G4Neutron::Neutron();
    nEl->AddDataSet(new G4NeutronElasticXS());
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(nEl, neutron);
    if( param->ApplyFactorXS() ) {
      nEl->MultiplyCrossSectionBy( param->XSFactorNucleonElastic() );
    }
  }
}
