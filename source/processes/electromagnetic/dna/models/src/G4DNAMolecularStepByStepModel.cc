#include "G4DNAMolecularStepByStepModel.hh"
#include "G4VDNAReactionModel.hh"

G4DNAMolecularStepByStepModel::G4DNAMolecularStepByStepModel(const G4String& name) :
    G4VITModel(name),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fReactionTable))
{
    fTimeStepper = new G4DNAMoleculeEncounterStepper();
    fReactionProcess = new G4DNAMolecularReaction();

    fType1 = G4Molecule::ITType();
    fType2 = G4Molecule::ITType();
    fReactionModel = 0;
}

G4DNAMolecularStepByStepModel::~G4DNAMolecularStepByStepModel()
{
    if(fReactionModel) delete fReactionModel;
}

G4DNAMolecularStepByStepModel::G4DNAMolecularStepByStepModel(const G4DNAMolecularStepByStepModel& right):
    G4VITModel(right),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fReactionTable))
{
        fReactionTable = right.fReactionTable;
        if(right.fReactionModel)
        {
            fReactionModel = right.fReactionModel->Clone();
            ((G4DNAMolecularReaction*)  fReactionProcess)->SetReactionModel(fReactionModel);
            ((G4DNAMoleculeEncounterStepper*) 	fTimeStepper)->SetReactionModel(fReactionModel);
        }
        else fReactionModel = 0;
}

void G4DNAMolecularStepByStepModel::Initialize()
{
    fReactionModel ->SetReactionTable((const G4DNAMolecularReactionTable*)fReactionTable);
    G4VITModel::Initialize();
}

void G4DNAMolecularStepByStepModel::PrintInfo()
{
    G4cout << "DNAMolecularStepByStepModel will be used" << G4endl;
}
