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
#include "G4EmDNAChemistry.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ProcessManager.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElectronSolvatation.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"

#include "G4DNAMolecularDissociation.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4VDNAReactionModel.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"

/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessTable.hh"
#include "G4DNASecondOrderReaction.hh"
#include "G4MolecularConfiguration.hh"
/****/

G4EmDNAChemistry::G4EmDNAChemistry() : G4VUserChemistryList()
{
}

G4EmDNAChemistry::~G4EmDNAChemistry()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry::ConstructMolecule()
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

	G4MoleculeTable::Instance()->CreateMoleculeModel("H3Op", G4H3O::Definition(), 1);
	G4Molecule* OHm = G4MoleculeTable::Instance()->CreateMoleculeModel("OHm", // just a tag to store and retrieve from G4MoleculeTable
			G4OH::Definition(), -1, 5.0e-9*(m2/s));
	OHm->SetMass(17.0079*g/Avogadro * c_squared);
	G4MoleculeTable::Instance()->CreateMoleculeModel("OH", G4OH::Definition());
	G4MoleculeTable::Instance()->CreateMoleculeModel("e_aq", G4Electron_aq::Definition());
	G4MoleculeTable::Instance()->CreateMoleculeModel("H", G4Hydrogen::Definition());
	G4MoleculeTable::Instance()->CreateMoleculeModel("H2", G4H2::Definition());
}

void G4EmDNAChemistry::ConstructDissociationChannels()
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
    G4MolecularDissociationChannel* decCh1;
    G4MolecularDissociationChannel* decCh2;

    G4ElectronOccupancy* occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

    //////////////////////////////////////////////////////////
    //            EXCITATIONS                               //
    //////////////////////////////////////////////////////////
    G4DNAWaterExcitationStructure waterExcitation;
    //--------------------------------------------------------
    //---------------Excitation on the fifth layer------------

    decCh1 = new G4MolecularDissociationChannel("A^1B_1_Relaxation");
    decCh2 = new G4MolecularDissociationChannel("A^1B_1_DissociativeDecay");
    //Decay 1 : OH + H
    decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
    decCh1->SetProbability(0.35);
    decCh1 -> SetDisplacementType(G4DNAWaterDissociationDisplacer::NoDisplacement);

    decCh2->AddProduct(OH);
    decCh2->AddProduct(H);
    decCh2->SetProbability(0.65);
    decCh2 -> SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

    water->AddExcitedState("A^1B_1");
    water->AddDecayChannel("A^1B_1",decCh1);
    water->AddDecayChannel("A^1B_1",decCh2);

    occ->RemoveElectron(4,1); // this is the transition form ground state to
    occ->AddElectron(5,1);  // the first unoccupied orbital: A^1B_1

    water->AddeConfToExcitedState("A^1B_1", *occ);

    //--------------------------------------------------------
    //---------------Excitation on the fourth layer-----------
    decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relaxation_Channel");
    decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociativeDecay");
    G4MolecularDissociationChannel* decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIonisation_Channel");

    water->AddExcitedState("B^1A_1");

    //Decay 1 : energy
    decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
    decCh1->SetProbability(0.3);

    //Decay 2 : 2OH + H_2
    decCh2->AddProduct(H2);
    decCh2->AddProduct(OH);
    decCh2->AddProduct(OH);
    decCh2->SetProbability(0.15);
    decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

    //Decay 3 : OH + H_3Op + e_aq
    decCh3->AddProduct(OH);
    decCh3->AddProduct(H3O);
    decCh3->AddProduct(e_aq);
    decCh3->SetProbability(0.55);
    decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

	 *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(3); // this is the transition form ground state to
    occ->AddElectron(5,1);  // the first unoccupied orbital: B^1A_1

    water->AddeConfToExcitedState("B^1A_1", *occ);
    water->AddDecayChannel("B^1A_1",decCh1);
    water->AddDecayChannel("B^1A_1",decCh2);
    water->AddDecayChannel("B^1A_1",decCh3);

    //-------------------------------------------------------
    //-------------------Excitation of 3rd layer-----------------
    decCh1 = new G4MolecularDissociationChannel("Excitation3rdLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDissociationChannel("Excitation3rdLayer_Relaxation_Channel");

    water->AddExcitedState("Excitation3rdLayer");

    //Decay channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);

    decCh1->SetProbability(0.5);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

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
    decCh1 = new G4MolecularDissociationChannel("Excitation2ndLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDissociationChannel("Excitation2ndLayer_Relaxation_Channel");

    water->AddExcitedState("Excitation2ndLayer");

    //Decay Channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);

    decCh1->SetProbability(0.5);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

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
    decCh1 = new G4MolecularDissociationChannel("Excitation1stLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDissociationChannel("Excitation1stLayer_Relaxation_Channel");

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
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

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

    decCh1 = new G4MolecularDissociationChannel("Ionisation_Channel");

    //Decay Channel 1 : : OH + H_3Op
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(OH);
    decCh1->SetProbability(1);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay);

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
    decCh1 = new G4MolecularDissociationChannel("DissociativeAttachment");

    //Decay 1 : 2OH + H_2
    decCh1->AddProduct(H2);
    decCh1->AddProduct(OHm);
    decCh1->AddProduct(OH);
    decCh1->SetProbability(1);
    //decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::DissociativeAttachment);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->AddElectron(5,1); // H_2O^-
    water->AddExcitedState("DissociativeAttachment");
    water->AddeConfToExcitedState("DissociativeAttachment", *occ);
    water->AddDecayChannel("DissociativeAttachment",decCh1);

    delete occ;
}

void G4EmDNAChemistry::ConstructReactionTable(G4DNAMolecularReactionTable* theReactionTable)
{
    //G4DNAMolecularReactionTable* theReactionTable = G4DNAMolecularReactionTable::GetReactionTable();

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
    // e_aq + e_aq + 2H2O -> H2 + 2OH-
    G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(0.5e10*(1e-3*m3/(mole*s)), e_aq,e_aq);
    reactionData -> AddProduct(OHm);
    reactionData -> AddProduct(OHm);
    reactionData -> AddProduct(H2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // e_aq + *OH -> OH-
    reactionData = new G4DNAMolecularReactionData(2.95e10*(1e-3*m3/(mole*s)), e_aq, OH);
    reactionData->AddProduct(OHm);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // e_aq + H* + H2O -> H2 + OH-
    reactionData = new G4DNAMolecularReactionData(2.65e10*(1e-3*m3/(mole*s)), e_aq, H);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(H2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // e_aq + H3O+ -> H* + H2O
    reactionData = new G4DNAMolecularReactionData(2.11e10*(1e-3*m3/(mole*s)), e_aq, H3Op);
    reactionData->AddProduct(H);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // e_aq + H2O2 -> OH- + *OH
    reactionData = new G4DNAMolecularReactionData(1.41e10*(1e-3*m3/(mole*s)), e_aq, H2O2);
    reactionData->AddProduct(OHm);
    reactionData->AddProduct(OH);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // *OH + *OH -> H2O2
    reactionData = new G4DNAMolecularReactionData(0.44e10*(1e-3*m3/(mole*s)), OH, OH);
    reactionData->AddProduct(H2O2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // *OH + *H -> H2O
    theReactionTable -> SetReaction(1.44e10*(1e-3*m3/(mole*s)), OH, H);
    //------------------------------------------------------------------
    // *H + *H -> H2
    reactionData = new G4DNAMolecularReactionData(1.20e10*(1e-3*m3/(mole*s)), H, H);
    reactionData->AddProduct(H2);
    theReactionTable -> SetReaction(reactionData);
    //------------------------------------------------------------------
    // H3O+ + OH- -> 2H2O
    theReactionTable -> SetReaction(1.43e11*(1e-3*m3/(mole*s)), H3Op, OHm);
    //------------------------------------------------------------------
}

void G4EmDNAChemistry::ConstructProcess()
{
	G4VProcess* process =
			G4ProcessTable::GetProcessTable()->
			FindProcess("e-_G4DNAVibExcitation",
			"e-");

        if(process)
	{
	 G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*) process;
	 G4DNASancheExcitationModel* sancheExcitationMod =
			new G4DNASancheExcitationModel();
	 vibExcitation->SetEmModel(sancheExcitationMod);
	 sancheExcitationMod->ExtendLowEnergyLimit(0.025*eV);
	}

	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	
	// *** Electron Solvatation ***
	ph->RegisterProcess(
			new G4DNAElectronSolvatation("e-_G4DNAElectronSolvatation"),
			G4Electron::Definition());

	G4MoleculeTable* theMoleculeTable = G4MoleculeTable::Instance();
	G4MoleculeDefinitionIterator iterator =
			theMoleculeTable ->GetDefintionIterator();
	iterator.reset();
	while(iterator())
	{
		G4MoleculeDefinition* moleculeDef = iterator.value();

		if(moleculeDef != G4H2O::Definition())
		{
		//	G4cout << "Brownian motion added for : "
		//	<< moleculeDef->GetName() << G4endl;
			G4DNABrownianTransportation* brown =
					new G4DNABrownianTransportation();
			//   brown->SetVerboseLevel(4);
			ph->RegisterProcess(brown, moleculeDef);
		}
		else
		{
			G4DNAMolecularDissociation* dissociationProcess =
					new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
			dissociationProcess -> SetDecayDisplacer(moleculeDef,
					new G4DNAWaterDissociationDisplacer);
			dissociationProcess -> SetVerboseLevel(1);
			ph->RegisterProcess(dissociationProcess, moleculeDef);
		}
		// Warning : end of particles and processes are needed by EM Physics builders
	}
}

void G4EmDNAChemistry::ConstructTimeStepModel(
		G4DNAMolecularReactionTable* reactionTable)
{

	//=========================================
	// Diffusion controlled reaction model
	//=========================================
	/**
	 * The reaction model defines how to compute the reaction range between molecules
	 */

	G4VDNAReactionModel* reactionRadiusComputer =
			new G4DNASmoluchowskiReactionModel();
	reactionTable->PrintTable(reactionRadiusComputer);

	/**
	 * The StepByStep model tells the step manager how to behave before and after each step,
	 * how to compute the time steps.
	 */

	G4DNAMolecularStepByStepModel* stepByStep =
			new G4DNAMolecularStepByStepModel();
	stepByStep->SetReactionModel(reactionRadiusComputer);
	//((G4DNAMoleculeEncounterStepper*) stepByStep->GetTimeStepper())->SetVerbose(5);

	RegisterTimeStepModel(stepByStep, 0);
}
