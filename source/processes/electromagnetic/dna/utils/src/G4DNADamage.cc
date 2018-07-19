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
// $Id: G4DNADamage.cc 85152 2014-10-24 16:07:37Z matkara $
//
#include "G4DNADamage.hh"
#include "G4UnitsTable.hh"

G4ThreadLocal G4DNADamage* G4DNADamage::fpInstance(0);

G4DNAIndirectHit::G4DNAIndirectHit(const G4String& baseName,
                                   const G4Molecule* molecule,
                                   const G4ThreeVector& position,
                                   G4double time) :
    G4VDNAHit(), fpMolecule(molecule)
{
  fBaseName = baseName;
  fPosition = position;
  fTime = time;
}

G4DNAIndirectHit::~G4DNAIndirectHit()
{
  if (fpMolecule) delete fpMolecule;
  fpMolecule = 0;
}

void G4DNAIndirectHit::Print()
{
  G4cout << "Reaction : " << fpMolecule->GetName() << " + " << fBaseName
         << " at position : " << G4BestUnit(fPosition, "Length")
         << " and time : " << G4BestUnit(fTime, "Time") << G4endl;
}

G4DNADamage* G4DNADamage::Instance()
{
  if (!fpInstance) new G4DNADamage();

  return fpInstance;
}

G4DNADamage::G4DNADamage()
{
  fJustCountDamage = false;
  fNIndirectDamage = 0;
  fpInstance = this;
}

G4DNADamage::~G4DNADamage()
{
  for (int i = 0; i < (int) fIndirectHits.size(); i++)
  {
    if (fIndirectHits[i]) delete fIndirectHits[i];
  }
  fIndirectHits.clear();
}

void G4DNADamage::DeleteInstance()
{
  if (fpInstance) delete fpInstance;
  fpInstance = 0;
}

void G4DNADamage::Reset()
{
  fNIndirectDamage = 0;
  for (int i = 0; i < (int) fIndirectHits.size(); i++)
  {
    if (fIndirectHits[i]) delete fIndirectHits[i];
  }
  fIndirectHits.clear();
}

void G4DNADamage::AddIndirectDamage(const G4String& baseName,
                                    const G4Molecule* molecule,
                                    const G4ThreeVector& position,
                                    G4double time)
{
  if (fJustCountDamage)
  {
    fNIndirectDamage++;
    return;
  }

  G4DNAIndirectHit* indirectHit = 0;
  std::map<G4Molecule, const G4Molecule*>::iterator it = 
   fMolMap.find(*molecule);

  if (it == fMolMap.end())
  {
    G4Molecule* mol(0);
    fMolMap[*molecule] = (mol = new G4Molecule(*molecule));
    indirectHit = new G4DNAIndirectHit(baseName, mol, position, time);
  }
  else
  {
    indirectHit = new G4DNAIndirectHit(baseName, it->second, position, time);
  }
  fIndirectHits.push_back(indirectHit);
}
