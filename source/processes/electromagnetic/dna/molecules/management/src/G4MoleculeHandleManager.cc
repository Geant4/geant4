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
// $Id: G4MoleculeHandleManager.cc 84858 2014-10-21 16:08:22Z gcosmo $
//
#include "G4MoleculeHandleManager.hh"
#include "G4Molecule.hh"

using namespace std;

G4ThreadLocal G4MoleculeHandleManager* G4MoleculeHandleManager::fInstance(0);

G4MoleculeHandleManager::G4MoleculeHandleManager()
{
  //    G4cout << "G4MoleculeHandleManager::G4MoleculeHandleManager()" << G4endl;
}

G4bool G4MoleculeHandleManager::CompMoleculePointer::operator()(const G4Molecule* mol1,
                                                                const G4Molecule* mol2) const
{
  return (*mol1) < (*mol2);
}

G4MoleculeHandleManager::~G4MoleculeHandleManager()
{
  if (fMoleculeHandle.empty() == false)
  {
    MoleculeHandleMap::iterator it = fMoleculeHandle.begin();
    for (; it != fMoleculeHandle.end(); it++)
    {
      it->second.reset();
    }
  }
}

void G4MoleculeHandleManager::DeleteInstance()
{
  if (fInstance)
  {
    delete fInstance;
    fInstance = 0;
  }
}

G4MoleculeHandleManager* G4MoleculeHandleManager::Instance()
{
  if (!fInstance)
  {
    fInstance = new G4MoleculeHandleManager;
  }
  return fInstance;
}

G4MoleculeHandle G4MoleculeHandleManager::GetMoleculeHandle(const G4Molecule* molecule)
{
  MoleculeHandleMap::iterator it = fMoleculeHandle.find(molecule);
  G4MoleculeHandle molHandle;

  if (it != fMoleculeHandle.end())
  {
    molHandle = G4MoleculeHandle(it->second);
  }
  else
  {
    molHandle = G4MoleculeHandle(molecule);
    fMoleculeHandle.insert(make_pair(molecule, G4MoleculeHandle(molHandle)));
  }

  return molHandle;
}
