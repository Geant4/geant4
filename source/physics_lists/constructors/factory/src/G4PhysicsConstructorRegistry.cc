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
// $Id: 
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4PhysicsConstructorRegistry
//
// Author  W. Pokorski  21.09.2012
//
// Modifications:
//

#include "G4ios.hh"
#include <iomanip>

#include "G4PhysicsConstructorRegistry.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4PhysicsConstructorFactory.hh"

G4ThreadLocal G4PhysicsConstructorRegistry* G4PhysicsConstructorRegistry::theInstance = 0;

G4PhysicsConstructorRegistry* G4PhysicsConstructorRegistry::Instance()
{
  if(0 == theInstance) {
    static G4ThreadLocal G4PhysicsConstructorRegistry *manager_G4MT_TLS_ = 0 ; if (!manager_G4MT_TLS_) manager_G4MT_TLS_ = new  G4PhysicsConstructorRegistry  ;  G4PhysicsConstructorRegistry &manager = *manager_G4MT_TLS_;
    theInstance = &manager;
  }
  return theInstance;
}

G4PhysicsConstructorRegistry::G4PhysicsConstructorRegistry()
{}

G4PhysicsConstructorRegistry::~G4PhysicsConstructorRegistry()
{
  Clean();
}

void G4PhysicsConstructorRegistry::Clean()
{
  size_t n = physConstr.size(); 
  if(n > 0) {
    for (size_t i=0; i<n; ++i) {
      if(physConstr[i]) {
	G4VPhysicsConstructor* p = physConstr[i];
	physConstr[i] = 0;
	delete p;
      }
    }
    physConstr.clear();
  }
}

void G4PhysicsConstructorRegistry::Register(G4VPhysicsConstructor* p)
{
  if(!p) return;
  size_t n = physConstr.size(); 
  if(n > 0) {
    for (size_t i=0; i<n; ++i) {
      if(physConstr[i] == p) { return; }
    }
  }
  physConstr.push_back(p);
}

void G4PhysicsConstructorRegistry::DeRegister(G4VPhysicsConstructor* p)
{
  if ( !p ) return;
  size_t n = physConstr.size(); 
  if ( n > 0 ) {
    for (size_t i=0; i<n; ++i) {
      if ( physConstr[i] == p ) {
        physConstr[i] = 0;
	return;
      }
    }
  }
}

void G4PhysicsConstructorRegistry::AddFactory(G4String name, G4VBasePhysConstrFactory* factory)
{
  factories[name] = factory;
}

G4VPhysicsConstructor* G4PhysicsConstructorRegistry::GetPhysicsConstructor(const G4String& name)
{
  // check if factory exists...
  //
  if (factories.find(name)!=factories.end())
    {
        // we could store the list of called factories in some vector and
        // before returning we can could first check if this physics constructor was already instantiated
        // if yes, we can throw an exception saying that this physics can been already registered
        
      return factories[name]->Instantiate();
    }
  else
    {
      G4ExceptionDescription ED;
      ED << "The factory for the physics constructor ["<< name << "] does not exist!" << G4endl;
      G4Exception("G4PhysicsConstructorRegistry::GetPhysicsConstructor", "PhysicsList001", FatalException, ED);
      return 0;
    }
}

G4bool G4PhysicsConstructorRegistry::IsKnownPhysicsConstructor(const G4String& name)
{
  return ( factories.find(name) != factories.end() );
}


std::vector<G4String> G4PhysicsConstructorRegistry::AvailablePhysicsConstructors() const
{
  std::vector<G4String> avail;
  std::map<G4String,G4VBasePhysConstrFactory*>::const_iterator itr;
  for ( itr = factories.begin(); itr != factories.end(); ++itr ) {
    avail.push_back(itr->first);
  }

  return avail;
}

void G4PhysicsConstructorRegistry::PrintAvailablePhysicsConstructors() const
{
  std::vector<G4String> avail = AvailablePhysicsConstructors();
  G4cout << "G4VPhysicsConstructors in G4PhysicsConstructorRegistry are:"
         << G4endl;
  if ( avail.empty() ) G4cout << "... no registered processes" << G4endl;
  else {
    size_t n = avail.size();
    for (size_t i=0; i<n; ++i ) {
      G4cout << " [" << std::setw(3) << i << "] "
             << " \"" << avail[i] << "\"" << G4endl;
    }
  }
}

//
// External reference to phy ctor factories for running with 'static' 
// libraries to pull the references of the declared factories into the
// same compilation unit as the registry itself.
// No harm having them in the non-static case.
//
G4_REFERENCE_PHYSCONSTR_FACTORY(G4ChargeExchangePhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4DecayPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAChemistry);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option1);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option2);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option3);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option4);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option5);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option6);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option7);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_stationary);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_stationary_option2);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_stationary_option4);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_stationary_option6);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmExtraPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmLivermorePhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmLivermorePolarizedPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmLowEPPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmPenelopePhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysicsGS);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysicsSS);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysicsWVI);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option1);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option2);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option3);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option4);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4GenericBiasingPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronDElasticPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronElasticPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsHP);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsLEND);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsXS);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronHElasticPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronInelasticQBBC);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_ATL);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_HP);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_TRV);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTF_BIC);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsINCLXX);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsNuBeam);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BERT);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BERT_HP);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC_AllHP);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC_HP);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_FTFP_BERT);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGS_BIC);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4HadronPhysicsShielding);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4ImportanceBiasing);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4IonBinaryCascadePhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4IonElasticPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4IonINCLXXPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4IonPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4IonQMDPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4NeutronCrossSectionXS);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4NeutronTrackingCut);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4OpticalPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4ParallelWorldPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4RadioactiveDecayPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4SpinDecayPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4StepLimiterPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4StoppingPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY(G4WeightWindowBiasing);

