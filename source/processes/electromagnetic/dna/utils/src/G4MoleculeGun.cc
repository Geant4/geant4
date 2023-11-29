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
/*
 * MoleculeGun.cc
 *
 *  Created on: 29 janv. 2014
 *      Author: kara
 */

#include "G4MoleculeGun.hh"
#include "G4MoleculeTable.hh"
#include "G4Molecule.hh"
#include "G4MoleculeGunMessenger.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include <cassert>
#include "Randomize.hh"
//#include "G4MIWorkspace.hh"
#include "G4MolecularConfiguration.hh"

//------------------------------------------------------------------------------

template<>
void TG4MoleculeShoot<G4Track>::ShootAtRandomPosition(G4MoleculeGun* gun)
{
  G4ThreeVector positionInLocalCoordinate;

  for(G4int i = 0; i < fNumber; ++i)
  {
    RandomPosInBox(*fBoxSize, positionInLocalCoordinate);
    gun->BuildAndPushTrack(fMoleculeName,
                 fPosition + positionInLocalCoordinate,
                 fTime);
  }
}

//------------------------------------------------------------------------------

template<>
void TG4MoleculeShoot<G4Track>::ShootAtFixedPosition(G4MoleculeGun* gun)
{
  for(G4int i = 0; i < fNumber; ++i)
  {
    gun->BuildAndPushTrack(fMoleculeName, fPosition, fTime);
  }
}

//------------------------------------------------------------------------------

template<>
void TG4MoleculeShoot<G4Track>::Shoot(G4MoleculeGun* gun)
{
  if(fBoxSize) ShootAtRandomPosition(gun);
  else ShootAtFixedPosition(gun);
}

//------------------------------------------------------------------------------

//template<>
//void TG4MoleculeShoot<G4ContinuousMedium>::Shoot(G4MoleculeGun*)
//{
//    G4MolecularConfiguration* conf = G4MoleculeTable::Instance()
//        ->GetConfiguration(fMoleculeName);
//    G4MIWorkspace::GetWorldWorkspace()->GetSpeciesInCM().Add(conf,
//                                                             fNumber);
//}

//------------------------------------------------------------------------------

G4MoleculeGun::G4MoleculeGun()
{
  fpMessenger = new G4MoleculeGunMessenger(this);
}

//------------------------------------------------------------------------------

G4MoleculeGun::~G4MoleculeGun()
{
  if (fpMessenger) delete fpMessenger;
}

//------------------------------------------------------------------------------

void G4MoleculeGun::DefineTracks()
{
  for (std::size_t i = 0; i < fShoots.size(); i++)
  {
    fShoots[i]->Shoot(this);
  }
}

//------------------------------------------------------------------------------

void G4MoleculeGun::AddMolecule(const G4String& name,
                                const G4ThreeVector& position,
                                G4double time)
{
  G4shared_ptr<G4MoleculeShoot> shoot(new TG4MoleculeShoot<G4Track>());
  shoot->fMoleculeName = name;
  shoot->fPosition = position;
  shoot->fTime = time;
  fShoots.push_back(shoot);
}

//------------------------------------------------------------------------------

void G4MoleculeGun::AddNMolecules(std::size_t n,
                                  const G4String& moleculeName,
                                  const G4ThreeVector& position,
                                  G4double time)
{
  G4shared_ptr<G4MoleculeShoot> shoot(new TG4MoleculeShoot<G4Track>());
  shoot->fNumber = (G4int)n;
  shoot->fMoleculeName = moleculeName;
  shoot->fPosition = position;
  shoot->fTime = time;
  fShoots.push_back(shoot);
}

//------------------------------------------------------------------------------

void
G4MoleculeGun::AddMoleculesRandomPositionInBox(std::size_t n,
                                               const G4String& moleculeName,
                                               const G4ThreeVector& boxCenter,
                                               const G4ThreeVector& boxSize,
                                               G4double time)
{
  G4shared_ptr<G4MoleculeShoot> shoot(new TG4MoleculeShoot<G4Track>());
  shoot->fNumber = (G4int)n;
  shoot->fMoleculeName = moleculeName;
  shoot->fPosition = boxCenter;
  shoot->fBoxSize = new G4ThreeVector(boxSize);
  shoot->fTime = time;
  fShoots.push_back(shoot);
}

//------------------------------------------------------------------------------

void G4MoleculeGun::BuildAndPushTrack(const G4String& name,
                                      const G4ThreeVector& position,
                                      G4double time)
{
  G4MolecularConfiguration* conf =
      G4MoleculeTable::Instance()->GetConfiguration(name);
  assert(conf != 0);
  G4Molecule* molecule = new G4Molecule(conf);

  PushTrack(molecule->BuildTrack(time, position));
}

//------------------------------------------------------------------------------

void G4MoleculeGun::GetNameAndNumber(G4MoleculeGun::NameNumber& output)
{
  for(std::size_t i = 0 ; i < fShoots.size() ; ++i)
  {
    output[fShoots[i]->fMoleculeName]+=fShoots[i]->fNumber;
  }
}

//------------------------------------------------------------------------------

void G4MoleculeShoot::RandomPosInBox(const G4ThreeVector& boxSize,
                                     G4ThreeVector& output)
{
  output[0] = boxSize.x() * G4UniformRand() - boxSize.x()/2;
  output[1] = boxSize.y() * G4UniformRand() - boxSize.y()/2;
  output[2] = boxSize.z() * G4UniformRand() - boxSize.z()/2;
}

//------------------------------------------------------------------------------

G4MoleculeShoot::G4MoleculeShoot() :
    G4enable_shared_from_this<G4MoleculeShoot>()
{
  fMoleculeName = "";
  fTime = 0;
  fNumber = 1;
  fBoxSize = 0;
}

//------------------------------------------------------------------------------

G4MoleculeShoot::~G4MoleculeShoot()
{
  if(fBoxSize) delete fBoxSize;
}

//------------------------------------------------------------------------------

void
G4MoleculeGun::AddMoleculeShoot(G4shared_ptr<G4MoleculeShoot> shoot)
{
  fShoots.push_back(shoot);
}

void G4MoleculeGun::AddMoleculeInCMRepresentation(std::size_t n,
                                                  const G4String& moleculeName,
                                                  G4double time)
{
  G4shared_ptr<G4MoleculeShoot> shoot(new TG4MoleculeShoot<G4ContinuousMedium>());
  shoot->fNumber = (G4int)n;
  shoot->fMoleculeName = moleculeName;
  shoot->fTime = time;
  fShoots.push_back(shoot);
}
