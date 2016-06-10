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
#include "G4LEPTSElasticModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LEPTSElasticModel::G4LEPTSElasticModel(const G4String& modelName) 
  : G4VLEPTSModel( modelName )
{
  theXSType = XSElastic;
} // constructor


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LEPTSElasticModel::~G4LEPTSElasticModel() {
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4LEPTSElasticModel::Initialise(const G4ParticleDefinition* aParticle, 
                          const G4DataVector&)
{
  Init();
  BuildPhysicsTable( *aParticle );

  fParticleChangeForGamma = GetParticleChangeForGamma();
  
  // static const G4double   proton_mass_c2 = 938.272013 * MeV;
  // static const G4double  neutron_mass_c2 = 939.56536 * MeV;
  // static const G4double h2o_mass_c2 = 8*neutron_mass_c2 + 10*(proton_mass_c2 + electron_mass_c2);
  // G4cout << "mme " << h2o_mass_c2/MeV << " " << H2o_mass_c2/MeV << G4endl;
   
  const G4MaterialTable * materialTable = G4Material::GetMaterialTable() ;
  std::vector<G4Material*>::const_iterator matite;
  for( matite = materialTable->begin(); matite != materialTable->end(); matite++ ) {
    const G4Material * aMaterial = (*matite);
    theMassTarget[aMaterial]  = theMolecularMass[aMaterial] / (6.02214179e+23/CLHEP::mole) *CLHEP::c_light * CLHEP::c_light;
    theMassProjectile[aMaterial] = CLHEP::electron_mass_c2;

    if( verboseLevel >= 1) G4cout << "Material: " << aMaterial->GetName() << " MolecularMass: " << theMolecularMass[aMaterial]/(CLHEP::g/CLHEP::mole) << " g/mole "
       << " MTarget: " << theMassTarget[aMaterial]/CLHEP::MeV << " MeV" << G4endl;
  }
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LEPTSElasticModel::CrossSectionPerVolume(const G4Material* mate,
                                         const G4ParticleDefinition* aParticle,
                                         G4double kineticEnergy,
                                         G4double,
                                         G4double)
{
  if( kineticEnergy < theLowestEnergyLimit ) return DBL_MAX;
  return 1./GetMeanFreePath( mate, aParticle, kineticEnergy );

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4LEPTSElasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple* mateCuts,
                                 const G4DynamicParticle* aDynamicParticle,
                                 G4double,
                                 G4double)
{
  G4double P0KinEn = aDynamicParticle->GetKineticEnergy();
  G4ThreeVector P0Dir = aDynamicParticle->GetMomentumDirection();

  if( P0KinEn < theLowestEnergyLimit ) {
    fParticleChangeForGamma->ProposeMomentumDirection( P0Dir );
    fParticleChangeForGamma->SetProposedKineticEnergy( 0.);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit( P0KinEn);
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    if( verboseLevel > 2 ) G4cout << " ENERGY LOW " << P0KinEn - theLowestEnergyLimit << G4endl;
    return;
  }

  //-  G4ParticleDefinition * particleDefDef = aTrack.GetDefinition();
  //-  G4String partName = particleDefDef->GetParticleName();

  //  G4ThreeVector pos, pos0, dpos;

  //-  G4StepPoG4int * PostPoG4int = aStep.GetPostStepPoG4int();
  //-  G4ThreeVector r = PostPoG4int->GetPosition();

  //TypeOfInteraction=-10;
  
  const G4Material* aMaterial = mateCuts->GetMaterial();
  G4double ang = SampleAngle(aMaterial, P0KinEn/CLHEP::eV, 0.0);
  G4ThreeVector P1Dir = SampleNewDirection(P0Dir, ang);
#ifdef DEBUG_LEPTS
  if( verboseLevel >= 2 ) G4cout << " G4LEPTSElasticModel::SampleSecondaries( P1Dir " << P1Dir << " P0Dir " << P0Dir << " ang " << ang << G4endl;
#endif

  //G4ThreeVector P1Dir = SampleNewDirection(P0Dir, P0KinEn/eV, 0.0);
  //G4double  Energylost1= ElasticEnergyTransferWater2(P0KinEn, ang);
  G4double  Energylost = EnergyTransfer(P0KinEn, ang, theMassTarget[aMaterial], theMassProjectile[aMaterial]);
  if( verboseLevel >= 3 ) G4cout << " ELASTIC Energylost "<< Energylost << " = " << P0KinEn << " " <<ang << " " << theMassTarget[aMaterial] << "  " << theMassProjectile[aMaterial] << G4endl; 

  G4double P1KinEn = P0KinEn - Energylost;
  if( verboseLevel >= 3 ) G4cout << " ELASTIC " << P1KinEn << " = " << P0KinEn << " - " << Energylost << G4endl;
#ifdef DEBUG_LEPTS
   if( verboseLevel >= 2 ) G4cout << " G4LEPTSElasticModel::SampleSecondaries( SetProposedKineticEnergy " << P1KinEn << " " << P0KinEn << " - " << Energylost << G4endl;
#endif
  fParticleChangeForGamma->ProposeMomentumDirection( P1Dir );
  fParticleChangeForGamma->SetProposedKineticEnergy( P1KinEn);
  fParticleChangeForGamma->ProposeLocalEnergyDeposit( Energylost);
  //G4cout << "elasticEnergyLost: " << Energylost << G4endl;

#ifdef DEBUG_LEPTS
   if( verboseLevel >= 2 ) G4cout << " G4LEPTSElasticModel::SampleSecondaries( ProposeMomentumDirection " <<  fParticleChangeForGamma->GetProposedMomentumDirection() << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LEPTSElasticModel::EnergyTransfer(G4double E, G4double ang, G4double MT, G4double MP) 
{
  G4double co = std::cos(ang);
  G4double si = std::sin(ang);

  G4double W = ( (E+MP)*si*si + MT - co*std::sqrt(MT*MT-MP*MP*si*si) ) * E*(E+2*MP)
    / ( std::pow((E+MP+MT),2) - E*co*co*(E+2*MP) );

  //G4double W2 = 2*MP/MT*(1-co)*E;
  //G4cout << "WWWWWWWWW: " << W/E << " " << E/W << " " << W2/W << G4endl;
  //G4cout << "Mm " << MT/MeV << " " << MP/MeV << G4endl;

  return W;
}

