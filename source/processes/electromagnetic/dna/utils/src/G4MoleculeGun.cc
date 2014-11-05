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
#include <cassert>

G4MoleculeGun::G4MoleculeGun()
{
  // TODO Auto-generated constructor stub

  fpMessenger = new G4MoleculeGunMessenger();
}

G4MoleculeGun::~G4MoleculeGun()
{
  // TODO Auto-generated destructor stub
  if (fpMessenger) delete fpMessenger;
}

void G4MoleculeGun::DefineTracks()
{
  fpMessenger->DefineTracks(this);

  for (size_t i = 0; i < fTracks.size(); i++)
  {
    PushTrack(fTracks[i]);
  }

  fTracks.clear();
}

void G4MoleculeGun::AddMolecule(const G4String& name,
                                const G4ThreeVector& position,
                                double time)
{
  G4Track* track = BuildTrack(name, position, time);
  fTracks.push_back(track);
}

void G4MoleculeGun::AddNMolecules(size_t n,
                                  const G4String& name,
                                  const G4ThreeVector& position,
                                  double time)
{
  for (size_t i = 0; i < n; i++)
  {
    AddMolecule(name, position, time);
  }
}

G4Track* G4MoleculeGun::BuildTrack(const G4String& name,
                                   const G4ThreeVector& position,
                                   double time)
{
  G4Molecule* model = G4MoleculeTable::Instance()->GetMoleculeModel(name);
  assert(model != 0);
  G4Molecule* molecule = new G4Molecule(*model);

  return molecule->BuildTrack(time, position);
}
