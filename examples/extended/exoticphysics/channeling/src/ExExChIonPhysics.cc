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

#include "ExExChIonPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"

using namespace std;

// factory
#include "G4PhysicsConstructorFactory.hh"
//

// Wrapper
#include "XWrapperDiscreteProcess.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(ExExChIonPhysics);

G4ThreadLocal G4bool ExExChIonPhysics::wasActivated = false;
G4ThreadLocal G4BinaryLightIonReaction* ExExChIonPhysics::theIonBC = 0;
G4ThreadLocal G4HadronicInteraction*  ExExChIonPhysics::theFTFP = 0;
G4ThreadLocal G4VCrossSectionDataSet*
    ExExChIonPhysics::theNuclNuclData = 0;
G4ThreadLocal G4VComponentCrossSection*
    ExExChIonPhysics::theGGNuclNuclXS = 0;
G4ThreadLocal G4FTFBuilder* ExExChIonPhysics::theBuilder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChIonPhysics::ExExChIonPhysics(G4int ver)
: G4VPhysicsConstructor("ionInelasticFTFP_BIC"),verbose(ver)
{
    SetPhysicsType(bIons);
    if(verbose > 1) { G4cout << "### ExExChIonPhysics" << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChIonPhysics::ExExChIonPhysics(const G4String& nname)
: G4VPhysicsConstructor(nname),verbose(1)
{
    SetPhysicsType(bIons);
    if(verbose > 1) { G4cout << "### ExExChIonPhysics" << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChIonPhysics::~ExExChIonPhysics()
{
    //Explictly setting pointers to zero is actually needed.
    //These are static variables, in case we restart threads
    //we need to re-create objects
    delete theBuilder; theBuilder = 0;
    delete theGGNuclNuclXS; theGGNuclNuclXS = 0;
    delete theNuclNuclData; theNuclNuclData = 0;
    delete theIonBC; theIonBC = 0;
    delete theFTFP; theFTFP = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChIonPhysics::ConstructParticle()
{
    //  Construct ions
    G4IonConstructor pConstructor;
    pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChIonPhysics::ConstructProcess()
{
    if(wasActivated) { return; }
    wasActivated = true;
    
    G4double emax = 100.*TeV;
    
    G4ExcitationHandler* handler = new G4ExcitationHandler();
    G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(handler);
    
    // Binary Cascade
    theIonBC = new G4BinaryLightIonReaction(thePreCompound);
    theIonBC->SetMinEnergy(0.0);
    theIonBC->SetMaxEnergy(4*GeV);
    
    // FTFP
    theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
    theFTFP = theBuilder->GetModel();
    theFTFP->SetMinEnergy(2*GeV);
    theFTFP->SetMaxEnergy(emax);
    
    theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc();
    theNuclNuclData = new G4CrossSectionInelastic(theGGNuclNuclXS);
    
    AddProcess("dInelastic", G4Deuteron::Deuteron(),false);
    AddProcess("tInelastic",G4Triton::Triton(),false);
    AddProcess("He3Inelastic",G4He3::He3(),true);
    AddProcess("alphaInelastic", G4Alpha::Alpha(),true);
    AddProcess("ionInelastic",G4GenericIon::GenericIon(),true);
    
    if(verbose > 1) {
        G4cout << "ExExChIonPhysics::ConstructProcess done! "
        << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChIonPhysics::AddProcess(const G4String& name,
                                        G4ParticleDefinition* part,
                                        G4bool )//isIon)
{
    G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess(name, part);
    G4ProcessManager* pManager = part->GetProcessManager();
    
    
    hadi->AddDataSet(theNuclNuclData);
    
    hadi->RegisterMe(theIonBC);
    hadi->RegisterMe(theFTFP);
    
    XWrapperDiscreteProcess* hadi_wrapper = new XWrapperDiscreteProcess();
    hadi_wrapper->RegisterProcess(hadi,1);
    pManager->AddDiscreteProcess(hadi_wrapper);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
