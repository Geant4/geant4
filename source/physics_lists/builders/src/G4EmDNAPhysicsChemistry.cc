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
#include "G4EmDNAPhysicsChemistry.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4DNAMolecularDecayDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ProcessManager.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAElectronSolvatation.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

#include "G4DNAMolecularDecay.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4VDNAReactionModel.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"

// Warning : the following is needed in order to use EM Physics builders
// e+
#include "G4Positron.hh"
#include "G4eMultipleScattering.hh"
#include "G4UrbanMscModel95.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
// gamma
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
// end of warning

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysicsChemistry);

using namespace std;

G4EmDNAPhysicsChemistry::G4EmDNAPhysicsChemistry(G4int ver)
    : G4VPhysicsConstructor("G4EmDNAPhysicsChemistry"), verbose(ver)
{
    SetPhysicsType(bElectromagnetic);
}

G4EmDNAPhysicsChemistry::~G4EmDNAPhysicsChemistry()
{
    /** It should not be useful to delete those singletons
      * since they are based on auto_ptr. But to make sure
      * that it won't cause any trouble, we explicitly call
      * the following methods responsible for the destruction
      * of the singletons.
      */
    G4DNAChemistryManager::DeleteInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysicsChemistry::ConstructParticle()
{
    // bosons
    G4Gamma::Gamma();

    // leptons
    G4Electron::Electron();
    G4Positron::Positron();

    // baryons
    G4Proton::Proton();

    G4GenericIon::GenericIonDefinition();

    G4DNAGenericIonsManager * genericIonsManager;
    genericIonsManager=G4DNAGenericIonsManager::Instance();
    genericIonsManager->GetIon("alpha++");
    genericIonsManager->GetIon("alpha+");
    genericIonsManager->GetIon("helium");
    genericIonsManager->GetIon("hydrogen");
}
void G4EmDNAPhysicsChemistry::ConstructMolecules()
{
    //-----------------------------------
    // Create the definition
    G4H2O::Definition();
    G4Hydrogen::Definition();
    G4H3O::Definition();
    G4OH::Definition();
    G4Electron_aq::Definition();
    G4H2O2::Definition();
    G4H2::Definition();

    ConstructDecayChannels() ;
}

void G4EmDNAPhysicsChemistry::ConstructDecayChannels()
{
    //-----------------------------------
    //Create the dynamic objects
    G4Molecule * H3O = new G4Molecule(G4H3O::Definition());
    H3O -> RemoveElectron(4,1);

    G4Molecule * OH = new G4Molecule(G4OH::Definition());
    G4Molecule * e_aq = new G4Molecule(G4Electron_aq::Definition());

    G4Molecule * H = new G4Molecule(G4Hydrogen::Definition());
    G4Molecule * H2 = new G4Molecule(G4H2::Definition());

    G4Molecule* OHm = new G4Molecule(G4OH::Definition());
    OHm->AddElectron(3);
    OHm->SetMass(17.0079*g/Avogadro * c_squared);
    OHm->SetDiffusionCoefficient(5.0e-9*(m2/s));

    //-------------------------------------
    //Define the decay channels
    G4MoleculeDefinition* water = G4H2O::Definition();
    G4MolecularDecayChannel* decCh1;
    G4MolecularDecayChannel* decCh2;

    G4ElectronOccupancy* occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

    //////////////////////////////////////////////////////////
    //            EXCITATIONS                               //
    //////////////////////////////////////////////////////////
    G4DNAWaterExcitationStructure waterExcitation;
    //--------------------------------------------------------
    //---------------Excitation on the fifth layer------------

    decCh1 = new G4MolecularDecayChannel("A^1B_1_Relaxation");
    decCh2 = new G4MolecularDecayChannel("A^1B_1_DissociativeDecay");
    //Decay 1 : OH + H
    decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
    decCh1->SetProbability(0.35);
    decCh1 -> SetDisplacementType(G4DNAMolecularDecayDisplacer::NoDisplacement);

    decCh2->AddProduct(OH);
    decCh2->AddProduct(H);
    decCh2->SetProbability(0.65);
    decCh2 -> SetDisplacementType(G4DNAMolecularDecayDisplacer::A1B1_DissociationDecay);

    water->AddExcitedState("A^1B_1");
    water->AddDecayChannel("A^1B_1",decCh1);
    water->AddDecayChannel("A^1B_1",decCh2);

    occ->RemoveElectron(4,1); // this is the transition form ground state to
    occ->AddElectron(5,1);  // the first unoccupied orbital: A^1B_1

    water->AddeConfToExcitedState("A^1B_1", *occ);

    //--------------------------------------------------------
    //---------------Excitation on the fourth layer-----------
    decCh1 = new G4MolecularDecayChannel("B^1A_1_Relaxation_Channel");
    decCh2 = new G4MolecularDecayChannel("B^1A_1_DissociativeDecay");
    G4MolecularDecayChannel* decCh3 = new G4MolecularDecayChannel("B^1A_1_AutoIonisation_Channel");

    water->AddExcitedState("B^1A_1");

    //Decay 1 : energy
    decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
    decCh1->SetProbability(0.3);

    //Decay 2 : 2OH + H_2
    decCh2->AddProduct(H2);
    decCh2->AddProduct(OH);
    decCh2->AddProduct(OH);
    decCh2->SetProbability(0.15);
    decCh2->SetDisplacementType(G4DNAMolecularDecayDisplacer::B1A1_DissociationDecay);

    //Decay 3 : OH + H_3Op + e_aq
    decCh3->AddProduct(OH);
    decCh3->AddProduct(H3O);
    decCh3->AddProduct(e_aq);
    decCh3->SetProbability(0.55);
    decCh3->SetDisplacementType(G4DNAMolecularDecayDisplacer::AutoIonisation);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(3); // this is the transition form ground state to
    occ->AddElectron(5,1);  // the first unoccupied orbital: B^1A_1

    water->AddeConfToExcitedState("B^1A_1", *occ);
    water->AddDecayChannel("B^1A_1",decCh1);
    water->AddDecayChannel("B^1A_1",decCh2);
    water->AddDecayChannel("B^1A_1",decCh3);

    //-------------------------------------------------------
    //-------------------Excitation of 3rd layer-----------------
    decCh1 = new G4MolecularDecayChannel("Excitation3rdLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDecayChannel("Excitation3rdLayer_Relaxation_Channel");

    water->AddExcitedState("Excitation3rdLayer");

    //Decay channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);

    decCh1->SetProbability(0.5);
    decCh1->SetDisplacementType(G4DNAMolecularDecayDisplacer::AutoIonisation);

    //Decay channel 2 : energy
    decCh2->SetEnergy(waterExcitation.ExcitationEnergy(2));
    decCh2->SetProbability(0.5);

    //Electronic configuration of this decay
    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(2,1);
    occ->AddElectron(5,1);

    //Configure the water molecule
    water->AddeConfToExcitedState("Excitation3rdLayer", *occ);
    water->AddDecayChannel("Excitation3rdLayer",decCh1);
    water->AddDecayChannel("Excitation3rdLayer",decCh2);

    //-------------------------------------------------------
    //-------------------Excitation of 2nd layer-----------------
    decCh1 = new G4MolecularDecayChannel("Excitation2ndLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDecayChannel("Excitation2ndLayer_Relaxation_Channel");

    water->AddExcitedState("Excitation2ndLayer");

    //Decay Channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);

    decCh1->SetProbability(0.5);
    decCh1->SetDisplacementType(G4DNAMolecularDecayDisplacer::AutoIonisation);

    //Decay channel 2 : energy
    decCh2->SetEnergy(waterExcitation.ExcitationEnergy(3));
    decCh2->SetProbability(0.5);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(1,1);
    occ->AddElectron(5,1);

    water->AddeConfToExcitedState("Excitation2ndLayer", *occ);
    water->AddDecayChannel("Excitation2ndLayer",decCh1);
    water->AddDecayChannel("Excitation2ndLayer",decCh2);

    //-------------------------------------------------------
    //-------------------Excitation of 1st layer-----------------
    decCh1 = new G4MolecularDecayChannel("Excitation1stLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDecayChannel("Excitation1stLayer_Relaxation_Channel");

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(0,1);
    occ->AddElectron(5,1);

    water->AddExcitedState("Excitation1stLayer");
    water->AddeConfToExcitedState("Excitation1stLayer", *occ);

    //Decay Channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);
    decCh1->SetProbability(0.5);
    decCh1->SetDisplacementType(G4DNAMolecularDecayDisplacer::AutoIonisation);

    //Decay channel 2 : energy
    decCh2->SetEnergy(waterExcitation.ExcitationEnergy(4));
    decCh2->SetProbability(0.5);
    water->AddDecayChannel("Excitation1stLayer",decCh1);
    water->AddDecayChannel("Excitation1stLayer",decCh2);

    /////////////////////////////////////////////////////////
    //                  IONISATION                         //
    /////////////////////////////////////////////////////////
    //--------------------------------------------------------
    //------------------- Ionisation -------------------------
    water->AddExcitedState("Ionisation");

    decCh1 = new G4MolecularDecayChannel("Ionisation_Channel");

    //Decay Channel 1 : : OH + H_3Op
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(OH);
    decCh1->SetProbability(1);
    decCh1->SetDisplacementType(G4DNAMolecularDecayDisplacer::Ionisation_DissociationDecay);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(4,1); // this is a ionized h2O with a hole in its last orbital
    water->AddeConfToExcitedState("Ionisation", *occ);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(3,1);
    water->AddeConfToExcitedState("Ionisation", *occ);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(2,1);
    water->AddeConfToExcitedState("Ionisation", *occ);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(1,1);
    water->AddeConfToExcitedState("Ionisation", *occ);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(0,1);
    water->AddeConfToExcitedState("Ionisation", *occ);
    water->AddDecayChannel("Ionisation",decCh1);
    // to this electronic configuration should be associated a decay time of 10e-15 s should the process do it on the dynamic object? the dyn object.

    //////////////////////////////////////////////////////////
    //            Dissociative Attachment                   //
    //////////////////////////////////////////////////////////
    decCh1 = new G4MolecularDecayChannel("DissociativeAttachment");

    //Decay 1 : 2OH + H_2
    decCh1->AddProduct(H2);
    decCh1->AddProduct(OHm);
    decCh1->AddProduct(OH);
    decCh1->SetProbability(1);
    //decCh1->SetDisplacementType(G4DNAMolecularDecayDisplacer::DissociativeAttachment);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->AddElectron(5,1); // H_2O^-
    water->AddExcitedState("DissociativeAttachment");
    water->AddeConfToExcitedState("DissociativeAttachment", *occ);
    water->AddDecayChannel("DissociativeAttachment",decCh1);

    delete occ;
}

void G4EmDNAPhysicsChemistry::ConstructReactionTable()
{
    G4DNAMolecularReactionTable* theReactionTable = G4DNAMolecularReactionTable::GetReactionTable();

    G4Molecule* OHm = new G4Molecule(G4OH::Definition());
    OHm->AddElectron(3);
    OHm->SetMass(17.0079*g/Avogadro * c_squared);
    OHm->SetDiffusionCoefficient(5.0e-9*(m2/s));
    G4Molecule* OH = new G4Molecule(G4OH::Definition());
    G4Molecule* e_aq = new G4Molecule(G4Electron_aq::Definition());
    G4Molecule* H = new G4Molecule(G4Hydrogen::Definition());
    G4Molecule* H3Op = new G4Molecule(G4H3O::Definition());
    G4Molecule* H2O2 = new G4Molecule(G4H2O2::Definition());
    H3Op->RemoveElectron(4,1);

    G4Molecule* H2 = new G4Molecule(G4H2::Definition());

    //------------------------------------------------------------------
    //e_aq + e_aq + 2H2O -> H2 + 2OH-
    G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(0.5e10*(1e-3*m3/(mole*s)), e_aq,e_aq);
    reactionData -> AddProduct(OHm);
    reactionData -> AddProduct(OHm);
    reactionData -> AddProduct(H2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //e_aq + *OH -> OH-
    reactionData = new G4DNAMolecularReactionData(2.95e10*(1e-3*m3/(mole*s)), e_aq, OH);
    reactionData->AddProduct(OHm);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //e_aq + H* + H2O -> H2 + OH-
    reactionData = new G4DNAMolecularReactionData(2.65e10*(1e-3*m3/(mole*s)), e_aq, H);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(H2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //e_aq + H3O+ -> H* + H2O
    reactionData = new G4DNAMolecularReactionData(2.11e10*(1e-3*m3/(mole*s)), e_aq, H3Op);
    reactionData->AddProduct(H);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //e_aq + H2O2 -> OH- + *OH
    reactionData = new G4DNAMolecularReactionData(1.41e10*(1e-3*m3/(mole*s)), e_aq, H2O2);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(OH);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //*OH + *OH -> H2O2
    reactionData = new G4DNAMolecularReactionData(0.44e10*(1e-3*m3/(mole*s)), OH, OH);
    reactionData->AddProduct(H2O2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //*OH + *H -> H2O
    theReactionTable -> SetReaction(1.44e10*(1e-3*m3/(mole*s)), OH, H);
    //------------------------------------------------------------------
    //*H + *H -> H2
    reactionData = new G4DNAMolecularReactionData(1.20e10*(1e-3*m3/(mole*s)), H, H);
    reactionData->AddProduct(H2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    //H3O+ + OH- -> 2H2O
    theReactionTable -> SetReaction(1.43e11*(1e-3*m3/(mole*s)), H3Op, OHm);
    //------------------------------------------------------------------
    //OH + H2 -> H + H2O
    reactionData = new G4DNAMolecularReactionData(4.17e7*(1e-3*m3/(mole*s)), OH, H2);
    reactionData->AddProduct(H);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
}

void G4EmDNAPhysicsChemistry::ConstructProcess()
{

    ConstructMolecules();
    // Is placed outside ConstructParticles to avoid
    // the transportation process to be added to the molecules

    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    theParticleIterator->reset();
    while( (*theParticleIterator)() )
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4String particleName = particle->GetParticleName();
        G4String particleType = particle->GetParticleType();
        G4ProcessManager* pmanager = particle->GetProcessManager();

        if (particleName == "e-") {

            // *** Elastic scattering (two alternative models available) ***

            G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
            G4DNAChampionElasticModel * championElasticModel = new G4DNAChampionElasticModel();
            championElasticModel->SetKillBelowThreshold(0.*eV);
            theDNAElasticProcess->SetModel(championElasticModel);

            // or alternative model
            //theDNAElasticProcess->SetModel(new G4DNAScreenedRutherfordElasticModel());

            ph->RegisterProcess(theDNAElasticProcess, particle);

            // *** Excitation ***
            ph->RegisterProcess(new G4DNAExcitation("e-_G4DNAExcitation"), particle);

            // *** Ionisation ***
            ph->RegisterProcess(new G4DNAIonisation("e-_G4DNAIonisation"), particle);

            // *** Vibrational excitation ***
            G4DNAVibExcitation* vibExcitation = new G4DNAVibExcitation("e-_G4DNAVibExcitation");
            G4DNASancheExcitationModel* sancheExcitationMod = new G4DNASancheExcitationModel;
            vibExcitation -> SetModel(sancheExcitationMod);
            sancheExcitationMod->ExtendLowEnergyLimit(0.025*eV);
            ph->RegisterProcess(vibExcitation, particle);

            // *** Attachment ***
            ph->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), particle);

            // *** Electron Solvatation ***
            ph->RegisterProcess(new G4DNAElectronSolvatation("e-_G4DNAElectronSolvatation"), particle);
        } else if ( particleName == "proton" ) {
            ph->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);

        } else if ( particleName == "hydrogen" ) {
            ph->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), particle);

        } else if ( particleName == "alpha" ) {
            ph->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), particle);

        } else if ( particleName == "alpha+" ) {
            ph->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), particle);

        } else if ( particleName == "helium" ) {
            ph->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), particle);

        }
            // Extension to HZE proposed by Z. Francis

        else if ( particleName == "carbon" ) {
          ph->RegisterProcess(new G4DNAIonisation("carbon_G4DNAIonisation"), particle);

        } else if ( particleName == "nitrogen" ) {
          ph->RegisterProcess(new G4DNAIonisation("nitrogen_G4DNAIonisation"), particle);

        } else if ( particleName == "oxygen" ) {
          ph->RegisterProcess(new G4DNAIonisation("oxygen_G4DNAIonisation"), particle);

        } else if ( particleName == "iron" ) {
          ph->RegisterProcess(new G4DNAIonisation("iron_G4DNAIonisation"), particle);

        }

        // Warning : the following particles and processes are needed by EM Physics builders
        // They are taken from the default Livermore Physics list
        // These particles are currently not handled by Geant4-DNA

        // e+

        else if (particleName == "e+") {

            // Identical to G4EmStandardPhysics_option3

            G4eMultipleScattering* msc = new G4eMultipleScattering();
            msc->AddEmModel(0, new G4UrbanMscModel95());
            msc->SetStepLimitType(fUseDistanceToBoundary);
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.2, 100*um);

            ph->RegisterProcess(msc, particle);
            ph->RegisterProcess(eIoni, particle);
            ph->RegisterProcess(new G4eBremsstrahlung(), particle);
            ph->RegisterProcess(new G4eplusAnnihilation(), particle);

        }else if (particleName == "gamma") {
            G4double LivermoreHighEnergyLimit = GeV;

            G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
            G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel =
                    new G4LivermorePhotoElectricModel();
            theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
            ph->RegisterProcess(thePhotoElectricEffect, particle);

            G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
            G4LivermoreComptonModel* theLivermoreComptonModel =
                    new G4LivermoreComptonModel();
            theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
            ph->RegisterProcess(theComptonScattering, particle);

            G4GammaConversion* theGammaConversion = new G4GammaConversion();
            G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel =
                    new G4LivermoreGammaConversionModel();
            theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
            ph->RegisterProcess(theGammaConversion, particle);

            G4RayleighScattering* theRayleigh = new G4RayleighScattering();
            G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
            theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            theRayleigh->AddEmModel(0, theRayleighModel);
            ph->RegisterProcess(theRayleigh, particle);

        }
        else if(particleType == "Molecule" && particleName != "H_{2}O")
        {
            if  (pmanager==0)
            {
                particle->SetProcessManager(new G4ProcessManager(particle));
            }

            G4DNABrownianTransportation* brown = new G4DNABrownianTransportation();
            ph->RegisterProcess(brown, particle);
        }
        else if (particleName == "H_{2}O")
        {
            if  (pmanager==0)
            {
                particle->SetProcessManager(new G4ProcessManager(particle));
            }

            G4DNAMolecularDecay* decayProcess = new G4DNAMolecularDecay("H2O_DNAMolecularDecay");
            decayProcess -> SetDecayDisplacer(particle, new G4DNAMolecularDecayDisplacer);
            decayProcess -> SetVerboseLevel(1);
            ph->RegisterProcess(decayProcess, particle);

        }


        // Warning : end of particles and processes are needed by EM Physics builders

    }

    // Deexcitation
    //
    G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
    de->SetFluo(true);

    // Chemistry
    ConstructReactionTable();

    /**
      * Tells to the chemistry manager whether the chemistry
      * needs to be activated.
      * WARNING : if you don't use the chemistry do not activate it
      * otherwise it might generate memory leaks with tracks created but
      * not destroyed.
      */
    G4DNAChemistryManager::Instance()->SetChemistryActivation(true);
}
