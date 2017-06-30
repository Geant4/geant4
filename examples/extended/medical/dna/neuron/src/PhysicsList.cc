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
// $Id$
//
/// \file PhysicsList.cc 
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "CommandLineParser.hh"
#include "G4EmParameters.hh"
// for discrete physics constructors!
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAChemistry.hh"
// for condensed physics constructors!
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
// for hadronic physics constructors!
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4StoppingPhysics.hh"  
#include "G4IonBinaryCascadePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"  

// G4.10.3 and later version
#include "G4EmDNAPhysicsActivator.hh"

using namespace G4DNAPARSER ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
: G4VModularPhysicsList()
{
  double currentDefaultCut ; 
  currentDefaultCut   = 1.*micrometer ; 
  SetDefaultCutValue(currentDefaultCut); 
  SetVerboseLevel(1);
  // fixe lower limit for cut
  G4ProductionCutsTable::GetProductionCutsTable()->
                         SetEnergyRange(100*eV, 1*GeV);
 
  // Options of combination Geant4-DNA processes (Physics/Chemistry) 
  // with Standard and Hadronic Physics: 

  // a) DNAphysics and Livermore physics inside and outside neuron 
  if(CommandLineParser::GetParser()->GetCommandIfActive("-dnaliv"))
  {
    G4cout<< "Livermore + DNAphysics is activated!"<<G4endl;
    RegisterConstructor("G4EmLivermorePhysics");

  // G4.10.2 and before version
  //  G4EmParameters::Instance()->AddDNA("BoundingSlice","Opt0");
  //  G4EmParameters::Instance()->AddDNA("Soma","Opt5");
  //  G4EmParameters::Instance()->AddDNA("Dendrites","Opt5");
  //  G4EmParameters::Instance()->AddDNA("Axon","Opt5");

  // G4.10.3 and later version!
  // 'G4EmParameters' works together with 'G4EmDNAPhysicsActivator'
    RegisterPhysics(new G4EmDNAPhysicsActivator());
    G4EmParameters::Instance()->AddDNA("Soma","Opt5");
    G4EmParameters::Instance()->AddDNA("Dendrites","Opt5");
    G4EmParameters::Instance()->AddDNA("Axon","Opt5");
  }
  
  //  b) Livermore + DNAPhysics + DNAChemistry 
  else if(CommandLineParser::GetParser()->GetCommandIfActive("-dnachemON"))
  {
    G4cout<< "'Livermore + DNAphysics + DNAChemistry' is activated!"<<G4endl;
    RegisterConstructor("G4EmLivermorePhysics");

  // G4.10.2 and before version
  //  G4EmParameters::Instance()->AddDNA("BoundingSlice","Opt0");
  //  RegisterConstructor("G4EmDNAChemistry");

  // G4.10.3 and later version
    RegisterPhysics(new G4EmDNAPhysicsActivator());
    RegisterPhysics(new G4EmDNAChemistry());
    G4EmParameters::Instance()->AddDNA("Soma","Opt5");
    G4EmParameters::Instance()->AddDNA("Dendrites","Opt5");
    G4EmParameters::Instance()->AddDNA("Axon","Opt5");
  }

  //  d) "QGSP_BIC_EMY" package from hadrontherapy advanced example
  else if(CommandLineParser::GetParser()->GetCommandIfActive("-dnahad"))
  {
    G4cout
      << "HadronPhysics + LivermorePhysics + DNAPhysics is activated!"<<G4endl;
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
   // RegisterConstructor("G4RadioactiveDecayPhysics");  
 
    // EM physics: Livermore
    RegisterConstructor("G4EmLivermorePhysics");  

    // EM physics: DNAPhysics inside BoundingSliceVolume
    // G4.10.2 and before version
    //G4EmParameters::Instance()->AddDNA("BoundingSlice","Opt0");

    // G4.10.3 and later version
    RegisterPhysics(new G4EmDNAPhysicsActivator());
    G4EmParameters::Instance()->AddDNA("Soma","Opt5");
    G4EmParameters::Instance()->AddDNA("Dendrites","Opt5");
    G4EmParameters::Instance()->AddDNA("Axon","Opt5");
  }

  //  Only G4EmStandardPhysics or G4EmLivermorePhysics in all volume
  else 
  {
    G4cout<< "Only LivermorePhysics is activated!"<<G4endl;
    RegisterConstructor("G4EmLivermorePhysics"); 
    //RegisterConstructor("G4EmStandardPhysics_option4"); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterConstructor(const G4String& name)
{
  RegisterPhysics(G4PhysicsConstructorRegistry::Instance()->
      GetPhysicsConstructor(name));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
