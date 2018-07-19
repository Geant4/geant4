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
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4PhysicsListHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
: G4VUserPhysicsList()
{
    SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
    ConstructBosons();
    ConstructLeptons();
    ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructBosons()
{
    // gamma
    G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructLeptons()
{
    // leptons
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//DNA
#include "G4DNAGenericIonsManager.hh"
//ENDDNA

void PhysicsList::ConstructBarions()
{
    // baryons
    G4Proton::ProtonDefinition();
    G4GenericIon::GenericIonDefinition();

    // Geant4 DNA new particles
    G4DNAGenericIonsManager * genericIonsManager;
    genericIonsManager=G4DNAGenericIonsManager::Instance();
    genericIonsManager->GetIon("alpha++");
    genericIonsManager->GetIon("alpha+");
    genericIonsManager->GetIon("helium");
    genericIonsManager->GetIon("hydrogen");
    genericIonsManager->GetIon("carbon");
    genericIonsManager->GetIon("nitrogen");
    genericIonsManager->GetIon("oxygen");
    genericIonsManager->GetIon("iron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructProcess()
{
    AddTransportation();
    ConstructEM();
    ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Geant4-DNA MODELS

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAPTBElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"
#include "G4DNAPTBExcitationModel.hh"

#include "G4DNAIonisation.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNARuddIonisationModel.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"
#include "G4DNAPTBIonisationModel.hh"

#include "G4DNAChargeDecrease.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"

#include "G4DNAChargeIncrease.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAMeltonAttachmentModel.hh"

#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"

//

#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4DNAVacuumModel.hh"
#include "G4VEmModel.hh"
#include "G4DNAModelInterface.hh"

//#include "G4ElectronCapture.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOOG4DNAPTB0OOooo........oooOO0OOooo....

void PhysicsList::ConstructEM()
{
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    auto myParticleIterator=GetParticleIterator();
    myParticleIterator->reset();
    while( (*myParticleIterator)() )
    {
        G4ParticleDefinition* particle = myParticleIterator->value();
        G4String particleName = particle->GetParticleName();

        if(particleName == "e-")
        {
            // **********************************************************
            // Instanciate models
            // **********************************************************

            G4DNAScreenedRutherfordElasticModel* e_modelDNARutherfordElastic = 
            new G4DNAScreenedRutherfordElasticModel();
            G4DNAEmfietzoglouIonisationModel* e_modelDNAEmfietzoglouIonisation =
            new G4DNAEmfietzoglouIonisationModel();
            G4DNAEmfietzoglouExcitationModel* e_modelDNAEmfietzoglouExcitation =
            new G4DNAEmfietzoglouExcitationModel();
            //            G4DNAMeltonAttachmentModel* e_modelDNAAttachement= 
            new G4DNAMeltonAttachmentModel();
            //            G4DNASancheExcitationModel* e_modelDNAVibra= 
            new G4DNASancheExcitationModel();

            G4DNAPTBElasticModel* modelDNAPTBElastic = 
            new G4DNAPTBElasticModel("THF/TMP/PY", particle);
            G4DNAPTBIonisationModel* modelDNAPTBIonisation = 
            new G4DNAPTBIonisationModel("THF/TMP/PY",particle);
            G4DNAPTBExcitationModel* modelDNAPTBExcitation = 
            new G4DNAPTBExcitationModel("THF/TMP/PY",particle);


            // Register models in the model interface

            //for elastics in three types of materials
            G4DNAModelInterface* e_elasticInteraction = 
            new G4DNAModelInterface("e-_elastic_interaction");
            e_elasticInteraction->
            RegisterModel(e_modelDNARutherfordElastic, particle);
            e_elasticInteraction->RegisterModel(modelDNAPTBElastic);
            e_elasticInteraction->RegisterModel(new G4DNAVacuumModel());


            //for ionization in three types of material

            G4DNAModelInterface* e_ionisationInteraction= 
            new G4DNAModelInterface("e-_ionisation_interaction");
            e_ionisationInteraction->
            RegisterModel(e_modelDNAEmfietzoglouIonisation,particle);
            e_ionisationInteraction->RegisterModel(modelDNAPTBIonisation);
            e_ionisationInteraction->RegisterModel(new G4DNAVacuumModel());

            //for excitation in three types of material
            G4DNAModelInterface* e_excitationInteraction= 
            new G4DNAModelInterface("e-_excitation_interaction");
            e_excitationInteraction->
            RegisterModel(e_modelDNAEmfietzoglouExcitation,particle);
            e_excitationInteraction->RegisterModel(modelDNAPTBExcitation);
            e_excitationInteraction->RegisterModel(new G4DNAVacuumModel());

            //Instanciate Processes
            // Elastic
            G4DNAElastic* e_DNAElasticProcess = 
            new G4DNAElastic("e-_G4DNAElastic");
            e_DNAElasticProcess->SetEmModel(e_elasticInteraction);
            // Ionisation
            G4DNAIonisation* e_DNAIonisationProcess = 
            new G4DNAIonisation("e-_G4DNAIonisation");
            e_DNAIonisationProcess->SetEmModel(e_ionisationInteraction);
            // Excitation
            G4DNAExcitation* e_DNAExcitationProcess = 
            new G4DNAExcitation("e-_G4DNAExcitation");
            e_DNAExcitationProcess->SetEmModel(e_excitationInteraction);
        /* G4DNAAttachment* e_DNAAttachementProcess= 
        new G4DNAAttachment("e-_G4DNAAttachement");
            e_DNAAttachementProcess->SetEmModel(e_attachementInteraction);
            G4DNAVibExcitation* e_DNAVibraProcess= 
            new G4DNAVibExcitation("e-_G4DNAVibraExci");
            e_DNAVibraProcess->SetEmModel(e_VibInteraction);*/


            // **********************************************************
            // Add previous process
            // **********************************************************
            
            // *** Elastic ***
            ph->RegisterProcess(e_DNAElasticProcess, particle);
            // *** Excitation ***
            ph->RegisterProcess(e_DNAExcitationProcess, particle);
            // *** Ionisation ***
            ph->RegisterProcess(e_DNAIonisationProcess, particle);


        } /*else if ( particleName == "proton" ) {

            G4DNAPTBIonisationModel* p_modelDNAPTBIonisation = 
            new G4DNAPTBIonisationModel("THF/TMP/PY", particle);
            G4DNARuddIonisationModel* p_modelDNARuddIonisation = 
            new G4DNARuddIonisationModel();
            //       G4DNABornIonisationModel* p_modelDNABornIonisation =
             new G4DNABornIonisationModel();
            //        G4DNABornExcitationModel* p_modelDNABornExcitation = 
            new G4DNABornExcitationModel();

            G4DNAModelInterface* p_ionisationInteraction= 
            new G4DNAModelInterface("p_ionisation_interaction");
            p_ionisationInteraction->RegisterModel(p_modelDNAPTBIonisation);
            p_ionisationInteraction->RegisterModel(new G4DNAVacuumModel());
            p_ionisationInteraction->RegisterModel(p_modelDNARuddIonisation, particle);

            G4DNAIonisation* p_DNAIonisationProcess = 
            new G4DNAIonisation("p_G4DNAIonisation");
            p_DNAIonisationProcess->SetEmModel(p_ionisationInteraction);
            ph->RegisterProcess(p_DNAIonisationProcess, particle);

        }*/
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructGeneral()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetCuts()
{
//Set default cut value to 1 nm for all particles

    SetDefaultCutValue(0.000001);
}


