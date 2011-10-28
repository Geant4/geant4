#include "G4DNASancheSolvatationModel.hh"
#include "G4WaterExcitationStructure.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4Molecule.hh"
#include "G4Electron_aq.hh"
#include "G4ITManager.hh"
#include "G4ITStepManager.hh"
#include "G4DNAChemistryManager.hh"

G4DNASancheSolvatationModel::G4DNASancheSolvatationModel(const G4ParticleDefinition*,
        const G4String& nam):
    G4VEmModel(nam),isInitialised(false)
{
    verboseLevel = 0 ;
    SetLowEnergyLimit(0.025*eV);
    G4WaterExcitationStructure exStructure ;
    SetHighEnergyLimit(exStructure.ExcitationEnergy(0));
}

//______________________________________________________________________
G4DNASancheSolvatationModel::~G4DNASancheSolvatationModel()
{}

//______________________________________________________________________
void G4DNASancheSolvatationModel::Initialise(const G4ParticleDefinition*,
        const G4DataVector& /*cuts*/)
{
    if (verboseLevel > 1)
        G4cout << "Calling G4SancheSolvatationModel::Initialise()" << G4endl;

    if(!isInitialised)
    {
        isInitialised = true;

        G4WaterExcitationStructure exStructure ;
        SetHighEnergyLimit(exStructure.ExcitationEnergy(0));

        if(pParticleChange)
            fParticleChangeForGamma = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
        else
            fParticleChangeForGamma = new G4ParticleChangeForGamma();
    }

}

//______________________________________________________________________
G4double G4DNASancheSolvatationModel::CrossSectionPerVolume(const G4Material* material,
        const G4ParticleDefinition* particleDefinition,
        G4double ekin,
        G4double,
        G4double)
{
    if (verboseLevel > 1)
        G4cout << "Calling CrossSectionPerVolume() of G4SancheSolvatationModel" << G4endl;

    if (particleDefinition != G4Electron::ElectronDefinition())
    {        
        __Exception_Origin__
        G4String exceptionCode ("G4DNASancheSolvatationModel001");
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "Attempting to calculate cross section for wrong particle";
        G4Exception(exceptionOrigin.data(),exceptionCode.data(),
                    FatalErrorInArgument,exceptionDescription);
        return 0;
    }

    if (material->GetName() == "G4_WATER")
    {
        if (ekin <= HighEnergyLimit())
        {
            return DBL_MAX;
        }
    }
    return 0. ;
}

//______________________________________________________________________
G4ThreeVector G4DNASancheSolvatationModel::radialDistributionOfProducts(G4double expectationValue) const
{
    G4double sigma = sqrt(1.57)/2*expectationValue;

    G4double XValueForfMax = sqrt(2.*sigma*sigma);
    G4double fMaxValue = sqrt(2./3.14) * 1./(sigma*sigma*sigma) *
                         (XValueForfMax*XValueForfMax)*
                         exp(-1./2. * (XValueForfMax*XValueForfMax)/(sigma*sigma));

    G4double R;

    do
    {
        G4double aRandomfValue = fMaxValue*G4UniformRand();

        G4double sign;
        if(G4UniformRand() > 0.5)
        {
            sign = +1.;
        }
        else
        {
            sign = -1;
        }

        R = expectationValue + sign*3.*sigma* G4UniformRand();
        G4double f = sqrt(2./3.14) * 1/pow(sigma, 3) * R*R * exp(-1./2. * R*R/(sigma*sigma));

        if(aRandomfValue < f)
        {
            break;
        }
    }
    while(1);

    G4double costheta = (2.*G4UniformRand()-1.);
    G4double theta = acos (costheta);
    G4double phi = 2.*pi*G4UniformRand();

    G4double xDirection = R*cos(phi)* sin(theta);
    G4double yDirection = R*sin(theta)*sin(phi);
    G4double zDirection = R*costheta;
    G4ThreeVector RandDirection = G4ThreeVector(xDirection, yDirection, zDirection);

    return RandDirection;
}

//______________________________________________________________________
void G4DNASancheSolvatationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
        const G4MaterialCutsCouple* /*couple*/,
        const G4DynamicParticle* particle,
        G4double,
        G4double)
{
    if (verboseLevel > 0)
        G4cout << "Calling SampleSecondaries() of G4SancheSolvatationModel" << G4endl;

    G4double k = particle->GetKineticEnergy();

    if (k <= HighEnergyLimit())
    {
        if(G4DNAChemistryManager::Instance()->IsChemistryActived())
        {
            G4double r_mean =
                (-0.003*pow(k/eV,6) + 0.0749*pow(k/eV,5) - 0.7197*pow(k/eV,4)
                 + 3.1384*pow(k/eV,3) - 5.6926*pow(k/eV,2) + 5.6237*k/eV - 0.7883)*nanometer;

            G4ThreeVector displacement = radialDistributionOfProducts (r_mean);
            //______________________________________________________________
            G4Molecule* e_aq = new G4Molecule(G4Electron_aq::Definition());
            const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
            G4Track * e_aqTrack = e_aq->BuildTrack(picosecond,theIncomingTrack->GetPosition()+displacement);
            e_aqTrack -> SetTrackStatus(fAlive);

            G4ITStepManager::Instance()->PushTrack(e_aqTrack);
            G4ITManager<G4Molecule>::Instance()->Push(e_aqTrack);
        }

        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    }
}
