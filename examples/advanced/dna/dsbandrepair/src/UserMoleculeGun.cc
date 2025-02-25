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
/// \file UserMoleculeGun.cc
/// \brief Implementation of the UserMoleculeGun class

#include "UserMoleculeGun.hh"
#include "UserMolecule.hh"

#include "G4MoleculeTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VITTrackHolder.hh"
#include "G4H2O.hh"
#include "G4Track.hh"

UserMoleculeShoot::UserMoleculeShoot() :
G4enable_shared_from_this<UserMoleculeShoot>()
{
    fMoleculeName = "";
    fTime = 0;
    fNumber = 1;
    fBoxSize = 0;
    fCopyNumber = -1;
    fStrand = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UserMoleculeShoot::~UserMoleculeShoot()
{
    if(fBoxSize) delete fBoxSize;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<>
void TUserMoleculeShoot<G4Track>::ShootAtRandomPosition(UserMoleculeGun* gun)
{
    G4ThreeVector positionInLocalCoordinate;

    for(int i = 0; i < fNumber; ++i)
    {
        RandomPosInBox(*fBoxSize, positionInLocalCoordinate);
        if(fStrand<=0)
        gun->BuildAndPushTrack(fMoleculeName,
        fPosition + positionInLocalCoordinate,
        fTime);
    else
        gun->BuildAndPushTrack(fMoleculeName,
        fPosition + positionInLocalCoordinate,
        fTime,fCopyNumber,fStrand);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<>
void TUserMoleculeShoot<G4Track>::ShootAtFixedPosition(UserMoleculeGun* gun)
{
    for(int i = 0; i < fNumber; ++i)
    {
        if(fStrand<=0)
        gun->BuildAndPushTrack(fMoleculeName, fPosition, fTime);
        else
        gun->BuildAndPushTrack(fMoleculeName,fPosition,fTime,fCopyNumber,fStrand);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<>
void TUserMoleculeShoot<G4Track>::MyShoot(UserMoleculeGun* gun)
{
    if(fBoxSize) ShootAtRandomPosition(gun);
    else ShootAtFixedPosition(gun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UserMoleculeGun::DefineTracks()
{
    for (size_t i = 0; i < fShoots.size(); i++)
    {
        fShoots[i]->MyShoot(this);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UserMoleculeGun::AddMolecule(const G4String& name,
    const G4ThreeVector& position,
    G4double time,
    G4int copyNumber,
    G4int strand)
{
    G4shared_ptr<UserMoleculeShoot> shoot(new TUserMoleculeShoot<G4Track>());
    shoot->fMoleculeName = name;
    shoot->fPosition = position;
    shoot->fTime = time;
    shoot->fCopyNumber = copyNumber;
    shoot->fStrand = strand;
    fShoots.push_back(shoot);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UserMoleculeGun::AddWaterMolecule(const G4ThreeVector& position,
    G4int trackId,
    ElectronicModification elecModif,
    G4int electronicLevel)
{
    UserMolecule * H2O = new UserMolecule (G4H2O::Definition() );
    switch (elecModif)
    {
        case eDissociativeAttachment:
        H2O -> AddElectron(5,1);
        break;
        case eExcitedMolecule :
        H2O -> ExciteMolecule(electronicLevel);
        break;
        case eIonizedMolecule :
        H2O -> IonizeMolecule(electronicLevel);
        break;
    }
    
    G4Track * H2OTrack = H2O->BuildTrack(1*picosecond, position);
    H2OTrack->SetParentID(trackId);
    H2OTrack->SetTrackStatus(fStopButAlive);
    H2OTrack->SetKineticEnergy(0.);
    G4VITTrackHolder::Instance()->Push(H2OTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UserMoleculeGun::BuildAndPushTrack(const G4String& name,
    const G4ThreeVector& position,
    G4double time)
{
    G4MolecularConfiguration* conf =
    G4MoleculeTable::Instance()->GetConfiguration(name);
    assert(conf != 0);
    UserMolecule* molecule = new UserMolecule(conf);
    PushTrack(molecule->BuildTrack(time, position));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UserMoleculeGun::BuildAndPushTrack(const G4String& name,
    const G4ThreeVector& position,
    G4double time,
    G4int copyNumber,
    G4int strand)
{
    G4MoleculeDefinition* def = G4MoleculeTable::Instance()->GetMoleculeDefinition(name);
    assert(def != 0);
    UserMolecule* molecule = new UserMolecule(def);
    molecule->SetCopyNumber(copyNumber);
    molecule->SetStrand(strand);

    PushTrack(molecule->BuildTrack(time, position));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
