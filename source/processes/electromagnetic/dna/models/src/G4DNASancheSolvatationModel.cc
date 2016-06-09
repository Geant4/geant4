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
// $Id: G4DNASancheSolvatationModel.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNASancheSolvatationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

G4DNASancheSolvatationModel::G4DNASancheSolvatationModel(const G4ParticleDefinition*,
                                                         const G4String& nam):
    G4VEmModel(nam),fIsInitialised(false)
{
    fVerboseLevel = 0 ;
    SetLowEnergyLimit(0.);
    G4DNAWaterExcitationStructure exStructure ;
    SetHighEnergyLimit(exStructure.ExcitationEnergy(0));
    fParticleChangeForGamma = 0;
  //  fNistWater  = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    fpWaterDensity = 0;
}

//______________________________________________________________________
G4DNASancheSolvatationModel::~G4DNASancheSolvatationModel()
{}

//______________________________________________________________________
void G4DNASancheSolvatationModel::Initialise(const G4ParticleDefinition* particleDefinition,
                                             const G4DataVector&)
{
#ifdef G4VERBOSE
    if (fVerboseLevel)
        G4cout << "Calling G4SancheSolvatationModel::Initialise()" << G4endl;
#endif
    if (particleDefinition != G4Electron::ElectronDefinition())
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "Attempting to calculate cross section for wrong particle";
        G4Exception("G4DNASancheSolvatationModel::CrossSectionPerVolume","G4DNASancheSolvatationModel001",
                    FatalErrorInArgument,exceptionDescription);
        return;
    }

    if(!fIsInitialised)
    {
        fIsInitialised = true;
        fParticleChangeForGamma = GetParticleChangeForGamma();
    }

    // Initialize water density pointer
    fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));
}

//______________________________________________________________________
G4double G4DNASancheSolvatationModel::CrossSectionPerVolume(const G4Material* material,
                                                            const G4ParticleDefinition*,
                                                            G4double ekin,
                                                            G4double,
                                                            G4double)
{
#ifdef G4VERBOSE
    if (fVerboseLevel > 1)
        G4cout << "Calling CrossSectionPerVolume() of G4SancheSolvatationModel" << G4endl;
#endif

    if(ekin > HighEnergyLimit())
    {
        return 0.0;
    }

    G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

    if(waterDensity!= 0.0)
   //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
    {
        if (ekin <= HighEnergyLimit())
        {
            return DBL_MAX;
        }
    }
    return 0. ;
}

//______________________________________________________________________
G4ThreeVector G4DNASancheSolvatationModel::RadialDistributionOfProducts(G4double expectationValue) const
{
    G4double sigma = std::sqrt(1.57)/2*expectationValue;

    G4double XValueForfMax = std::sqrt(2.*sigma*sigma);
    G4double fMaxValue = std::sqrt(2./3.14) * 1./(sigma*sigma*sigma) *
            (XValueForfMax*XValueForfMax)*
            std::exp(-1./2. * (XValueForfMax*XValueForfMax)/(sigma*sigma));

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
        G4double f = std::sqrt(2./3.14) * 1/std::pow(sigma, 3) * R*R * std::exp(-1./2. * R*R/(sigma*sigma));

        if(aRandomfValue < f)
        {
            break;
        }
    }
    while(1);

    G4double costheta = (2.*G4UniformRand()-1.);
    G4double theta = std::acos (costheta);
    G4double phi = 2.*pi*G4UniformRand();

    G4double xDirection = R*std::cos(phi)* std::sin(theta);
    G4double yDirection = R*std::sin(theta)*std::sin(phi);
    G4double zDirection = R*costheta;
    G4ThreeVector RandDirection = G4ThreeVector(xDirection, yDirection, zDirection);

    return RandDirection;
}

//______________________________________________________________________
void G4DNASancheSolvatationModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                                    const G4MaterialCutsCouple*,
                                                    const G4DynamicParticle* particle,
                                                    G4double,
                                                    G4double)
{
#ifdef G4VERBOSE
    if (fVerboseLevel)
        G4cout << "Calling SampleSecondaries() of G4SancheSolvatationModel" << G4endl;
#endif
    G4double k = particle->GetKineticEnergy();

    if (k <= HighEnergyLimit())
    {
        G4double r_mean =
                (-0.003*std::pow(k/eV,6) + 0.0749*std::pow(k/eV,5) - 0.7197*std::pow(k/eV,4)
                 + 3.1384*std::pow(k/eV,3) - 5.6926*std::pow(k/eV,2) + 5.6237*k/eV - 0.7883)*nanometer;

        G4ThreeVector displacement = RadialDistributionOfProducts (r_mean);

        //______________________________________________________________
        const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4ThreeVector finalPosition(theIncomingTrack->GetPosition()+displacement);

        G4DNAChemistryManager::Instance()->CreateSolvatedElectron(theIncomingTrack,&finalPosition);

        fParticleChangeForGamma->SetProposedKineticEnergy(25.e-3*eV);
        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    }
}
