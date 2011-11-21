#include "G4DNATransformElectronModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"

G4DNATransformElectronModel::G4DNATransformElectronModel(const G4ParticleDefinition*,
                                                         const G4String& nam):
    G4VEmModel(nam),fIsInitialised(false)
{
    fVerboseLevel = 0 ;
    SetLowEnergyLimit(0.*eV);
    SetHighEnergyLimit(0.025*eV);
    fParticleChangeForGamma = 0;
}

//______________________________________________________________________
G4DNATransformElectronModel::~G4DNATransformElectronModel()
{}

//______________________________________________________________________
void G4DNATransformElectronModel::Initialise(const G4ParticleDefinition* particleDefinition,
                                             const G4DataVector&)
{
#ifdef G4VERBOSE
    if (fVerboseLevel)
        G4cout << "Calling G4DNATransformElectronModel::Initialise()" << G4endl;
#endif

    if (particleDefinition != G4Electron::ElectronDefinition())
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "Attempting to calculate cross section for wrong particle";
        G4Exception("G4DNATransformElectronModel::CrossSectionPerVolume","G4DNATransformElectronModel001",
                    FatalErrorInArgument,exceptionDescription);
        return;
    }

    if(!fIsInitialised)
    {
        fIsInitialised = true;
        fParticleChangeForGamma = GetParticleChangeForGamma();
        fNistWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    }
}

//______________________________________________________________________
G4double G4DNATransformElectronModel::CrossSectionPerVolume(const G4Material* material,
                                                            const G4ParticleDefinition*,
                                                            G4double ekin,
                                                            G4double,
                                                            G4double)
{
#if G4VERBOSE
    if (fVerboseLevel > 1)
        G4cout << "Calling CrossSectionPerVolume() of G4DNATransformElectronModel" << G4endl;
#endif

    if(ekin > HighEnergyLimit())
    {
        return 0.0;
    }

    if (material == fNistWater || material->GetBaseMaterial() == fNistWater)
    {
        if (ekin <= HighEnergyLimit())
        {
            return DBL_MAX;
        }
    }

    return 0.0 ;
}

//______________________________________________________________________
void G4DNATransformElectronModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                    const G4MaterialCutsCouple*,
                                                    const G4DynamicParticle* particle,
                                                    G4double,
                                                    G4double)
{
#if G4VERBOSE
    if (fVerboseLevel)
        G4cout << "Calling SampleSecondaries() of G4DNATransformElectronModel" << G4endl;
#endif

    G4double k = particle->GetKineticEnergy();

    if (k <= HighEnergyLimit())
    {
        const G4Track * track = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateSolvatedElectron(track);
        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    }
}
