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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "CommandLineParser.hh"
#include "G4EmParameters.hh"
// for discrete physics constructors!
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4EmDNAChemistry.hh"
#include "G4EmDNAChemistry_option1.hh"
// for condensed physics constructors!
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
// for hadronic physics constructors!
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4StoppingPhysics.hh"  
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"  

#include "G4EmDNAPhysicsActivator.hh"

using namespace G4DNAPARSER ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
  : G4VModularPhysicsList(),
    fEmPhysicsList(nullptr),
    fDNAActivator(nullptr),
    fEmDNAChemistryList(nullptr),
    fEmDNAChemistryList1(nullptr),
    fEmName(""), fHadronic(false)
{
  G4double currentDefaultCut = 1.*micrometer; 
  SetDefaultCutValue(currentDefaultCut); 
  SetVerboseLevel(1);
  // fixe lower limit for cut
  G4ProductionCutsTable::GetProductionCutsTable()->
                         SetEnergyRange(100*eV, 1*GeV);
 
  // Options of combination Geant4-DNA processes (Physics/Chemistry) 
  // with Standard and Hadronic Physics: 

  // a) DNAphysics and Livermore physics inside and outside neuron 
  G4cout<< "Livermore + DNAphysics is activated!"<<G4endl;
  RegisterConstructor("emlivermore");

  // 'G4EmParameters' works together with 'G4EmDNAPhysicsActivator'
  // VI: in this example Livermore is a default in any way
  //     DNA option4 is the default configuration as well
  fDNAActivator = new G4EmDNAPhysicsActivator();
  /*
  G4EmParameters::Instance()->AddDNA("Soma","Opt4");
  G4EmParameters::Instance()->AddDNA("Dendrites","Opt4");
  G4EmParameters::Instance()->AddDNA("Axon","Opt4");
  */
  //  b) Livermore + DNAPhysics + DNAChemistry 
  if(CommandLineParser::GetParser()->GetCommandIfActive("-dnachemON"))
  {
    G4cout<< "DNAChemistry is activated!"<<G4endl;
    RegisterPhysics(new G4EmDNAChemistry());
  }

  //  d) "QGSP_BIC_EMY" package from hadrontherapy advanced example
  if(CommandLineParser::GetParser()->GetCommandIfActive("-dnahad"))
  {
    G4cout << "QGSP_BIC is activated!"<<G4endl;
    RegisterConstructor("QGSP_BIC");
    fHadronic = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(fEmPhysicsList)       { fEmPhysicsList->ConstructParticle(); }
  if(fDNAActivator)        { fDNAActivator->ConstructParticle(); }
  if(fEmDNAChemistryList)  { fEmDNAChemistryList->ConstructParticle(); }
  if(fEmDNAChemistryList1) { fEmDNAChemistryList1->ConstructParticle(); }
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  if(fEmPhysicsList)       { fEmPhysicsList->ConstructProcess(); }
  if(fDNAActivator)        { fDNAActivator->ConstructProcess(); }
  if(fEmDNAChemistryList)  { fEmDNAChemistryList->ConstructProcess(); }
  if(fEmDNAChemistryList1) { fEmDNAChemistryList1->ConstructProcess(); }
  G4VModularPhysicsList::ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterConstructor(const G4String& name)
{
  if(name == fEmName) { return; }
  if (name == "emstandard_opt0"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();

  } else if (name == "emstandard_opt3"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "emstandard_opt4"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
    
  } else if (name == "empenelope"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();

  } else if (name == "emlivermore"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();

  } else if(name == "G4EmDNAPhysics") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option1") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option1(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option2") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option2(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option3") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option3(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option4") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option4(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option5") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option5(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option6") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option6(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option7") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option7(verboseLevel);
    fEmName = name;

  } else if(name == "G4EmDNAPhysics_option8") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option8(verboseLevel);
    fEmName = name;

  } else if(name == "QGSP_BIC") {
    if(fHadronic) { return; }
    // Hadron Elastic Physics
    RegisterConstructor("G4HadronElasticPhysics");
    // Hadron Inelastic Physics
    RegisterConstructor("G4HadronPhysicsQGSP_BIC"); 
    // Stopping
    RegisterConstructor("G4StoppingPhysics");
    // Ion Physics
    RegisterConstructor("G4IonBinaryCascadePhysics");
    // Gamma-Lepto nuclear 
    RegisterConstructor("G4EmExtraPhysics");
    // Limiters
    RegisterConstructor("G4NeutronTrackingCut"); 
    // Decay
    RegisterConstructor("G4DecayPhysics");
    // Radioactive decay
    RegisterConstructor("G4RadioactiveDecayPhysics");  

  } else if(name == "G4EmDNAChemistry") {
    if(fEmDNAChemistryList || fEmDNAChemistryList1) { return; }
    fEmDNAChemistryList = new G4EmDNAChemistry();

  } else if(name == "G4EmDNAChemistry_option1") {
    if(fEmDNAChemistryList || fEmDNAChemistryList1) { return; }
    fEmDNAChemistryList1 = new G4EmDNAChemistry_option1();

  } else {
    G4cout << "PhysicsList::RegisterConstructor: <" << name << ">"
           << " fails - name is not defined"
           << G4endl;    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
