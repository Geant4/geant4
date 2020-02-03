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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4CrossSectionDataSetRegistry
//
// Author  V.Ivanchenko  24.01.2009
//
// Modifications:
//

#include "G4ios.hh"

#include "G4CrossSectionDataSetRegistry.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4VComponentCrossSection.hh"
#include "G4CrossSectionFactory.hh"
#include "G4CrossSectionFactoryRegistry.hh"

// Neeed for running with 'static' libraries to pull the references of the 
// declared factories
G4_REFERENCE_XS_FACTORY(G4ChipsKaonMinusInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsKaonMinusElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsKaonPlusInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsKaonPlusElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsKaonZeroInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsKaonZeroElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsHyperonInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsHyperonElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsProtonInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsProtonElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsNeutronInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsNeutronElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsPionPlusInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsPionPlusElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsPionMinusInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsPionMinusElasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsAntiBaryonInelasticXS);
G4_REFERENCE_XS_FACTORY(G4ChipsAntiBaryonElasticXS);
G4_REFERENCE_XS_FACTORY(G4NucleonNuclearCrossSection);
G4_REFERENCE_XS_FACTORY(G4ElectroNuclearCrossSection);
G4_REFERENCE_XS_FACTORY(G4PhotoNuclearCrossSection);
G4_REFERENCE_XS_FACTORY(G4PiNuclearCrossSection);
G4_REFERENCE_XS_FACTORY(G4NeutronInelasticXS);
G4_REFERENCE_XS_FACTORY(G4NeutronElasticXS);
G4_REFERENCE_XS_FACTORY(G4NeutronCaptureXS);

G4ThreadLocal G4CrossSectionDataSetRegistry* G4CrossSectionDataSetRegistry::instance = nullptr;

G4CrossSectionDataSetRegistry* G4CrossSectionDataSetRegistry::Instance()
{
  if(!instance) {
    static G4ThreadLocalSingleton<G4CrossSectionDataSetRegistry> inst;
    instance = inst.Instance();
  }
  return instance;
}

G4CrossSectionDataSetRegistry::G4CrossSectionDataSetRegistry()
{}

G4CrossSectionDataSetRegistry::~G4CrossSectionDataSetRegistry()
{
  Clean();
}

void G4CrossSectionDataSetRegistry::Clean()
{
  size_t n = xSections.size(); 
  for (size_t i=0; i<n; ++i) {
    G4VCrossSectionDataSet* p = xSections[i];
    if(p) {
      //std::cout << "Clean x-section #" << i << std::endl;
      //std::cout << "     " << p->GetName() << std::endl;
      delete p;
      xSections[i] = nullptr;
    }
  }
  n = xComponents.size(); 
  for (size_t i=0; i<n; ++i) {
    G4VComponentCrossSection* p = xComponents[i];
    if(p) {
      //std::cout << "Clean component #" << i << std::endl;
      //std::cout << "     " << p->GetName() << std::endl;
      delete p;
      xComponents[i] = nullptr;
    }
  }
}

void G4CrossSectionDataSetRegistry::Register(G4VCrossSectionDataSet* p)
{
  if(!p) { return; }
  for (auto xsec : xSections) {
    if(xsec == p) { return; }
  }
  //G4cout << "Register x-section #" << xSections.size() 
  //<< "  " << p->GetName() << G4endl;
  xSections.push_back(p);
}

void G4CrossSectionDataSetRegistry::DeRegister(G4VCrossSectionDataSet* p)
{
  if(!p) { return; }
  size_t n = xSections.size(); 
  for (size_t i=0; i<n; ++i) {
    if(xSections[i] == p) {
      //std::cout << "DeRegister x-section #" << i << std::endl;
      //std:: cout << "        " << p->GetName() << std::endl;
      xSections[i] = nullptr;
      return;
    }
  }
}

void G4CrossSectionDataSetRegistry::Register(G4VComponentCrossSection* p)
{
  if(!p) { return; }
  for (auto xsec : xComponents) {
    if(xsec == p) { return; }
  }
  //G4cout << "Register component #" << xComponents.size() 
  //<< "  " << p->GetName() << G4endl;
  xComponents.push_back(p);
}

void G4CrossSectionDataSetRegistry::DeRegister(G4VComponentCrossSection* p)
{
  if(!p) { return; }
  size_t n = xComponents.size(); 
  for (size_t i=0; i<n; ++i) {
    if(xComponents[i] == p) {
      //std::cout << "DeRegister component #" << i << std::endl;
      //std::cout << "         " << p->GetName() << std::endl;
      xComponents[i] = nullptr;
      return;
    }
  }
}

void 
G4CrossSectionDataSetRegistry::DeleteComponent(G4VComponentCrossSection* p)
{
  if(!p) { return; }
  size_t n = xComponents.size(); 
  for (size_t i=0; i<n; ++i) {
    if(xComponents[i] == p) {
      //std::cout << "Delete component #" << i << std::endl;
      //std::cout << "      " << p->GetName() << std::endl;
      delete p;
      return;
    }
  }
}

G4VCrossSectionDataSet* 
G4CrossSectionDataSetRegistry::GetCrossSectionDataSet(const G4String& name, 
                                                      G4bool warning)
{
  for (auto xsec : xSections) {
    if(xsec && xsec->GetName() == name) { return xsec; }
  }
  // check if factory exists...
  //
  G4CrossSectionFactoryRegistry* factories = 
    G4CrossSectionFactoryRegistry::Instance();
  // This throws if factory is not found, add second parameter 
  // to false to avoid this
  G4VBaseXSFactory* factory = factories->GetFactory(name, warning );
  if ( factory ) {
    return factory->Instantiate();
  } else {
    G4VCrossSectionDataSet* ptr = nullptr;
    return ptr;
  }
}

G4VComponentCrossSection*
G4CrossSectionDataSetRegistry::GetComponentCrossSection(const G4String& name)
{
  for (auto xsec : xComponents) {
    if(xsec && xsec->GetName() == name) { return xsec; }
  }
  G4VComponentCrossSection* ptr = nullptr;
  return ptr;
}
