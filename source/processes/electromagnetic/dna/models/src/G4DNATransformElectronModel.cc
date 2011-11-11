#include "G4DNATransformElectronModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4Molecule.hh"
#include "G4Electron_aq.hh"
#include "G4ITStepManager.hh"
#include "G4ITManager.hh"
#include "G4DNAChemistryManager.hh"

G4DNATransformElectronModel::G4DNATransformElectronModel(const G4ParticleDefinition*,
        const G4String& nam):
    G4VEmModel(nam),isInitialised(false)
{
    verboseLevel = 0 ;
    SetLowEnergyLimit(0.*eV);
    SetHighEnergyLimit(0.025*eV);
    fParticleChangeForGamma = 0;
}

//______________________________________________________________________
G4DNATransformElectronModel::~G4DNATransformElectronModel()
{}

//______________________________________________________________________
void G4DNATransformElectronModel::Initialise(const G4ParticleDefinition* particleDefinition,
        const G4DataVector& /*cuts*/)
{
    if (verboseLevel > 1)
        G4cout << "Calling G4DNATransformElectronModel::Initialise()" << G4endl;

    if (particleDefinition != G4Electron::ElectronDefinition())
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "Attempting to calculate cross section for wrong particle";
        G4Exception("G4DNATransformElectronModel::CrossSectionPerVolume","G4DNATransformElectronModel001",
                    FatalErrorInArgument,exceptionDescription);
        return;
    }

    if(!isInitialised)
    {
        isInitialised = true;
        highEnergyLimit = 0.025*eV;
        SetHighEnergyLimit(highEnergyLimit);
        fParticleChangeForGamma = GetParticleChangeForGamma();
    }
}

//______________________________________________________________________
G4double G4DNATransformElectronModel::CrossSectionPerVolume(const G4Material* material,
        const G4ParticleDefinition* /*particleDefinition*/,
        G4double ekin,
        G4double,
        G4double)
{
    if (verboseLevel > 1)
        G4cout << "Calling CrossSectionPerVolume() of G4DNATransformElectronModel" << G4endl;

    if(ekin > highEnergyLimit) return 0;

    if (material->GetName() == "G4_WATER")
    {
        if (ekin <= highEnergyLimit)
        {
            return DBL_MAX;
        }
    }

    return -1 ;
}

//______________________________________________________________________
void G4DNATransformElectronModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
        const G4MaterialCutsCouple* /*couple*/,
        const G4DynamicParticle* particle,
        G4double,
        G4double)
{
    if (verboseLevel > 0)
        G4cout << "Calling SampleSecondaries() of G4DNATransformElectronModel" << G4endl;

    G4double k = particle->GetKineticEnergy();

    if (k <= highEnergyLimit)
    {
        const G4Track * track = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateSolvatedElectron(track);
        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    }
}
