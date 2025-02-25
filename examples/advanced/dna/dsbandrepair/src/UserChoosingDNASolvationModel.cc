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
//
/// \file UserChoosingDNASolvationModel.cc
/// \brief Implementation of the UserChoosingDNASolvationModel class

#include "UserChoosingDNASolvationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4TransportationManager.hh"
#include "G4ITNavigator.hh"
#include "G4Navigator.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename MODEL>
UserTDNAOneStepThermalizationModel<MODEL>::
UserTDNAOneStepThermalizationModel(const G4ParticleDefinition*,const G4String& nam) :
G4TDNAOneStepThermalizationModel<MODEL>(0,nam)
{
    G4cout << "Calling SampleSecondaries() of UserTDNAOneStepThermalizationModel for Solvation process!!! \n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename MODEL>
void UserTDNAOneStepThermalizationModel<MODEL>::
SampleSecondaries(std::vector<G4DynamicParticle*>*,
                  const G4MaterialCutsCouple*,
                  const G4DynamicParticle* particle,
                  G4double tmin,
                  G4double maxEnergy)
{
    G4TDNAOneStepThermalizationModel<MODEL>::SampleSecondaries(0,0,particle,tmin,maxEnergy);
    
    if (G4DNAChemistryManager::IsActivated()) return; // 
    G4double k = particle->GetKineticEnergy();
    G4double highE = G4TDNAOneStepThermalizationModel<MODEL>::HighEnergyLimit();
    if (k <= highE){
        G4ThreeVector displacement(0,0,0);
        G4TDNAOneStepThermalizationModel<MODEL>::GetPenetration(k, displacement);

        auto _ParticleChangeForGamma=G4TDNAOneStepThermalizationModel<MODEL>::GetParticleChangeForGamma();
        const G4Track * theIncomingTrack =_ParticleChangeForGamma->GetCurrentTrack();

        G4ThreeVector finalPosition(theIncomingTrack->GetPosition()+displacement);
        this->fpNavigator->SetWorldVolume(theIncomingTrack->GetTouchable()->GetVolume(theIncomingTrack->GetTouchable()->GetHistoryDepth()));

        G4double displacementMag = displacement.mag();
        G4double safety = DBL_MAX;
        G4ThreeVector direction = displacement/displacementMag;

        G4double mag_displacement = displacement.mag();
        G4ThreeVector displacement_direction = displacement/mag_displacement;


        this->fpNavigator->ResetHierarchyAndLocate(theIncomingTrack->GetPosition(),
                                            direction,
                                            *((G4TouchableHistory*)
                                            theIncomingTrack->GetTouchable()));
        this->fpNavigator->ComputeStep(theIncomingTrack->GetPosition(),
                                displacement/displacementMag,
                                displacementMag,
                                safety);

        if(safety <= displacementMag){
            finalPosition = theIncomingTrack->GetPosition()+ (displacement/displacementMag)*safety*0.80;
        }
        G4DNAChemistryManager::Instance()->CreateSolvatedElectron(theIncomingTrack,&finalPosition);
        _ParticleChangeForGamma->SetProposedKineticEnergy(25.e-3*eV);
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel* UserChoosingDNASolvationModel::UserCreate(const G4String& penetrationModel)
{
    G4String modelNamePrefix("DNAOneStepThermalizationModel_");
    
    if(penetrationModel == "Terrisol1990")
    {
        return new UserTDNAOneStepThermalizationModel<DNA::Penetration::Terrisol1990>
        (G4Electron::Definition(), modelNamePrefix + penetrationModel);
    }
    else if(penetrationModel == "Meesungnoen2002")
    {
        return new UserTDNAOneStepThermalizationModel<DNA::Penetration::Meesungnoen2002>
        (G4Electron::Definition(), modelNamePrefix + penetrationModel);
    }
    else if(penetrationModel == "Meesungnoen2002_amorphous")
    {
        return new UserTDNAOneStepThermalizationModel<DNA::Penetration::Meesungnoen2002_amorphous>
        (G4Electron::Definition(), modelNamePrefix + penetrationModel);
    }
    else if(penetrationModel == "Kreipl2009")
    {
        return new UserTDNAOneStepThermalizationModel<DNA::Penetration::Kreipl2009>
        (G4Electron::Definition(), modelNamePrefix + penetrationModel);
    }
    else if(penetrationModel == "Ritchie1994")
    {
        return new UserTDNAOneStepThermalizationModel<DNA::Penetration::Ritchie1994>
        (G4Electron::Definition(), modelNamePrefix + penetrationModel);
    }
    else
    {
        G4ExceptionDescription description;
        description << penetrationModel + " is not a valid model name.";
        G4Exception("UserChoosingDNASolvationModel::UserCreate",
                    "INVALID_ARGUMENT",
                    FatalErrorInArgument,
                    description,
                    "Options are: Terrisol1990, Meesungnoen2002, Ritchie1994.");
    }
    return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel* UserChoosingDNASolvationModel::UserGetMacroDefinedModel()
{
    auto dnaSubType = G4EmParameters::Instance()->DNAeSolvationSubType();
    
    switch(dnaSubType)
    {
    case fRitchie1994eSolvation:
        return UserCreate("Ritchie1994");
    case fTerrisol1990eSolvation:
        return UserCreate("Terrisol1990");
    case fKreipl2009eSolvation:
        return UserCreate("Kreipl2009");
    case fMeesungnoensolid2002eSolvation:
        return UserCreate("Meesungnoen2002_amorphous");
    case fMeesungnoen2002eSolvation:
    case fDNAUnknownModel:
        return UserCreate("Meesungnoen2002");
    default:
        G4ExceptionDescription msg;
        msg<<"The solvation parameter stored in G4EmParameters is unknown. "
            <<"Supported types are: fRitchie1994eSolvation, fTerrisol1990eSolvation, "
            <<"fMeesungnoen2002eSolvation.";
        G4Exception("UserChoosingDNASolvationModel::UserGetMacroDefinedModel",
                    "DnaSubType",
                    FatalErrorInArgument,msg);
    }

    return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......