#include "G4LEPTSAttachmentModel.hh"

// constructor
G4LEPTSAttachmentModel::G4LEPTSAttachmentModel(const G4String& modelName) 
  : G4VLEPTSModel( modelName )
{
  theXSType = XSAttachment;

} // constructor


G4LEPTSAttachmentModel::~G4LEPTSAttachmentModel() {
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4LEPTSAttachmentModel::Initialise(const G4ParticleDefinition* aParticle, 
                          const G4DataVector&)
{
  Init();
  BuildPhysicsTable( *aParticle );

  fParticleChangeForGamma = GetParticleChangeForGamma();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LEPTSAttachmentModel::CrossSectionPerVolume(const G4Material* mate,
                                         const G4ParticleDefinition* aParticle,
                                         G4double kineticEnergy,
                                         G4double,
                                         G4double)
{
  return 1./GetMeanFreePath( mate, aParticle, kineticEnergy );

}


void G4LEPTSAttachmentModel::SampleSecondaries(std::vector<G4DynamicParticle*>* ,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle* aDynamicParticle,
                                 G4double,
                                 G4double)
{
  G4double P0KinEn = aDynamicParticle->GetKineticEnergy();

  fParticleChangeForGamma->SetProposedKineticEnergy(0.);
  fParticleChangeForGamma->ProposeLocalEnergyDeposit (P0KinEn);
  fParticleChangeForGamma->ProposeTrackStatus( fStopAndKill);
 
}
