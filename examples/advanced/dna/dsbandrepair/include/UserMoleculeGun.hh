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
/// \file UserMoleculeGun.hh
/// \brief Definition of the UserMoleculeGun class

#ifndef UserMoleculeGun_h
#define UserMoleculeGun_h 1

#include "G4DNAChemistryManager.hh"
#include "G4MoleculeGun.hh"

class UserMoleculeGun;

class UserMoleculeShoot : public G4enable_shared_from_this<UserMoleculeShoot>, 
public G4MoleculeShoot
{
public:

    UserMoleculeShoot();
    ~UserMoleculeShoot() override;
    void Shoot(G4MoleculeGun*) override {};
    virtual void MyShoot(UserMoleculeGun*) = 0;

    template<typename TYPE> G4shared_ptr<UserMoleculeShoot> ChangeType();

    G4int fCopyNumber{-1};
    G4int fStrand{-1};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename TYPE>
class TUserMoleculeShoot : public UserMoleculeShoot
{
public:
    TUserMoleculeShoot() : UserMoleculeShoot(){;}
    ~TUserMoleculeShoot() override {;}
    void MyShoot(UserMoleculeGun*) override;
protected:
    void ShootAtRandomPosition(UserMoleculeGun*);
    void ShootAtFixedPosition(UserMoleculeGun*);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename TYPE>
G4shared_ptr<UserMoleculeShoot> UserMoleculeShoot::ChangeType()
{
    G4shared_ptr<UserMoleculeShoot> output(new TUserMoleculeShoot<TYPE>);
    output->fMoleculeName = fMoleculeName;
    output->fPosition = fPosition;
    output->fTime = fTime;
    output->fNumber = fNumber;
    output->fBoxSize = fBoxSize;
    output->fCopyNumber = fCopyNumber;
    output->fStrand = fStrand;
    return output;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UserMoleculeGun : public G4ITGun
{
public:
    UserMoleculeGun() = default;
    ~UserMoleculeGun() override = default;

    void DefineTracks() override;

    void AddMolecule(const G4String& name,
    const G4ThreeVector& position,
    G4double time,
    G4int copyNumber,
    G4int strand);

    void AddWaterMolecule(const G4ThreeVector& position,
        G4int trackId,
        ElectronicModification elecModif,
        G4int electronicLevel);

protected:
    void BuildAndPushTrack(const G4String& name,
    const G4ThreeVector& position,
    G4double time = 0);

    void BuildAndPushTrack(const G4String& name,
                            const G4ThreeVector& position,
                            G4double time,
                            G4int copyNumber,
                            G4int strand);

    std::vector<G4shared_ptr<UserMoleculeShoot> > fShoots;
friend class UserMoleculeShoot;
template<class T> friend class TUserMoleculeShoot;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
