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
/// file: ChromosomeMapper.cc
/// brief: Implementation of class to track chromosomes
#include "ChromosomeMapper.hh"

#include "VirtualChromosome.hh"
#include "ChromosomeMessenger.hh"

#include <fstream>
#include <sstream>
#include "DNAHashing.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeMapper::ChromosomeMapper()
  : fChromosomes({})
{
  fpChromosomeMessenger = new ChromosomeMessenger(this);
  fpChromosomeFactory   = new ChromosomeFactory();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeMapper::~ChromosomeMapper()
{
  delete fpChromosomeMessenger;
  delete fpChromosomeFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

[[maybe_unused]] G4String ChromosomeMapper::GetCurrentChromosomeKey(
  const G4ThreeVector& pos) const
{
  G4String key = "";
  for(const auto& Chromosome : fChromosomes)
  {
    if(Chromosome.second->PointInChromosome(pos))
    {
      //key = Chromosome.first;
      key = Chromosome.second->GetName();
      break;
    }
  }
  return key;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

[[maybe_unused]] VirtualChromosome* ChromosomeMapper::GetChromosome(
  const G4ThreeVector& pos) const
{
  for(const auto& fChromosome : fChromosomes)
  {
    if(fChromosome.second->PointInChromosome(pos))
    {
      return fChromosome.second;
    }
  }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

[[maybe_unused]] VirtualChromosome* ChromosomeMapper::GetChromosome(
  const G4String& key) const
{
  uint32_t key_i = G4::hashing::crc32::Hash(key);

  try
  {
    return fChromosomes.at(key_i);
  } catch(const std::out_of_range& oor)
  {
    G4cout << "Chromosome does not exist with key: " << key << G4endl;
    return nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeMapper::AddChromosome(const G4String& key,
                                     const std::vector<G4String>& commands)
{
  uint32_t key_i = G4::hashing::crc32::Hash(key);

  auto it = fChromosomes.find(key_i);

  if(it == fChromosomes.end())
  {
    auto* newChromosome =
      ChromosomeFactory::MakeChromosome(key, commands);
    fChromosomes.emplace(key_i, newChromosome);
  }
  else
  {
    G4ExceptionDescription errmsg;
    errmsg << "ChromosomeMapper:: "
           << "Chromosome already exists with key: " << key << G4endl;
    G4Exception("ChromosomeMapper::AddChromosome", "", FatalException, errmsg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeMapper::SavePlotData(const G4String& filename)
{
  std::fstream fs(filename,
                  std::fstream::in | std::fstream::out | std::fstream::trunc);
  for(auto& fChromosome : fChromosomes)
  {
    fs << fChromosome.second->Plot();
  }
  fs.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

[[maybe_unused]] std::vector<G4String> ChromosomeMapper::GetChromosomeKeys()
                                                            const
{
  std::vector<G4String> keys;
  for(const auto& fChromosome : fChromosomes)
  {
    keys.push_back(fChromosome.second->GetName());
  }
  return keys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeMapper::Test()
{ fpChromosomeFactory->Test();}
