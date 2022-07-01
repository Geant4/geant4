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
/// \file ChromosomeHit.cc
/// \brief Hit class for a hit interacting with a DNA molecule

#include <utility>

#include "ChromosomeHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<ChromosomeHit>* MolecularChromosomeHitAllocator =
  nullptr;

ChromosomeHit::ChromosomeHit(G4String  name)
  : G4VHit()
  , fName(std::move(name))
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeHit::~ChromosomeHit() = default;

ChromosomeHit::ChromosomeHit(const ChromosomeHit& right)
  : G4VHit()
  , fName(right.GetName())
  , fEdepChromosome(right.GetChromosomeEdep())
  , fEdepDNA(right.GetDNAEdep())
{
  // consider setters?
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const ChromosomeHit& ChromosomeHit::operator=(const ChromosomeHit& right)
{
  auto* newHit = new ChromosomeHit(right);
  return *newHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ChromosomeHit::operator==(const ChromosomeHit& right) const
{
  return (this == &right) ? 1 : 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
