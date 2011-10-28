#include "G4DNADiffusionControlledModel.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4DNAMolecularReactionTable.hh"

pthread_mutex_t G4DNADiffusionControlledModel::fMutex = PTHREAD_MUTEX_INITIALIZER;

G4DNADiffusionControlledModel::G4DNADiffusionControlledModel() : G4VDNAReactionModel()
{
    fReactionData = 0 ;
}

G4DNADiffusionControlledModel::G4DNADiffusionControlledModel(const G4DNADiffusionControlledModel& __right) :
    G4VDNAReactionModel(__right)
{
//    G4cout << "G4DiffusionControlledModel::G4DiffusionControlledModel(const G4DiffusionControlledModel& right)"
//        << G4endl;
    fReactionData = 0 ;
}
G4DNADiffusionControlledModel::~G4DNADiffusionControlledModel()
{
    fReactionData = 0 ;
}

void G4DNADiffusionControlledModel::Initialise(const G4Molecule* __molecule, const G4Track&)
{
    fReactionData = fReactionTable->GetReactionData(__molecule);
}

void G4DNADiffusionControlledModel::InitialiseToPrint(const G4Molecule* __molecule)
{
    fReactionData = fReactionTable->GetReactionData(__molecule);
}

G4double G4DNADiffusionControlledModel::GetReactionRadius(const G4Molecule* mol1,
                                                               const G4Molecule* mol2)
{
    //pthread_mutex_lock(&fMutex);
    G4double __output = fReactionTable -> GetReactionData(mol1,mol2)->GetReducedReactionRadius();
    //pthread_mutex_unlock(&fMutex);

    return  __output ;
}

G4double G4DNADiffusionControlledModel::GetReactionRadius(const G4int __i)
{
    //pthread_mutex_lock(&fMutex);
    G4double __output = (*fReactionData)[__i] -> GetReducedReactionRadius();
    //pthread_mutex_unlock(&fMutex);

    return  __output ;
}

G4bool G4DNADiffusionControlledModel::FindReaction(const G4Track& __trackA,
        const G4Track& __trackB,
        const G4double __R,
        G4double& __r,
        const G4bool __alongStepReaction)
{
    G4double __postStepSeparation = 0;
    bool __do_break = false ;
    G4double __R2 = __R*__R ;
    int k = 0 ;

    for(; k < 3 ; k++)
    {
        __postStepSeparation += pow(__trackA.GetPosition()[k] - __trackB.GetPosition()[k],2);

        if(__postStepSeparation > __R2)
        {
            __do_break = true  ;
            break ;
        }
    }

    if(__do_break == false)
    {
         // The loop was not break
         // => __r^2 < __R^2
        __r = sqrt(__postStepSeparation);
        return true;
    }
    else if(__alongStepReaction == true)
    {
        //G4cout << "alongStepReaction==true" << G4endl;
        //Along step cheack and
        // the loop has break

        // Continue loop
        for(; k < 3 ; k++)
        {
            __postStepSeparation += pow(__trackA.GetPosition()[k] - __trackB.GetPosition()[k],2);
        }
        // Use Green approach : the Brownian bridge
        __r = (__postStepSeparation = sqrt(__postStepSeparation) );

        G4Molecule* __moleculeA = GetMolecule(__trackA);
        G4Molecule* __moleculeB = GetMolecule(__trackB);

        G4double __D = __moleculeA->GetDiffusionCoefficient() + __moleculeB->GetDiffusionCoefficient();

        G4ThreeVector __preStepPositionA = __trackA.GetStep()->GetPreStepPoint() ->GetPosition();
        G4ThreeVector __preStepPositionB = __trackB.GetStep()->GetPreStepPoint() ->GetPosition();

        if(__preStepPositionA == __trackA.GetPosition())
        G4Exception();

        G4double __preStepSeparation = (__preStepPositionA - __preStepPositionB).mag();

        G4double __probabiltyOfEncounter =  exp(-(__preStepSeparation - __R)*(__postStepSeparation - __R)
                                        / (__D* (__trackB.GetStep()->GetDeltaTime())));
        G4double __selectedPOE = G4UniformRand();

        if(__selectedPOE<=__probabiltyOfEncounter)  return true;
    }

    return false ;
}
