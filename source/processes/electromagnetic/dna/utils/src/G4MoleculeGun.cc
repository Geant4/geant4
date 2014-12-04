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
