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
// G4FermiBreakUp alternative de-excitation model
// by A. Novikov (January 2025)
//

#include "G4FermiIntegerPartition.hh"

using namespace fbu;

G4FermiIntegerPartition::G4FermiIntegerPartition(uint32_t number, uint32_t termsCount,
                                                 uint32_t base)
  : number_(number), termsCount_(termsCount), base_(base)
{}

G4FermiIntegerPartition::Iterator G4FermiIntegerPartition::begin() const
{
  return {number_, termsCount_, base_};
}

G4FermiIntegerPartition::Iterator G4FermiIntegerPartition::end() const
{
  return {};
}

/////////////////////////////////// ITERATOR //////////////////////////////

G4FermiIntegerPartition::Iterator::pointer G4FermiIntegerPartition::Iterator::operator->() const
{
  return &partition_;
}

G4FermiIntegerPartition::Iterator::reference G4FermiIntegerPartition::Iterator::operator*() const
{
  return partition_;
}

G4FermiIntegerPartition::Iterator& G4FermiIntegerPartition::Iterator::operator++()
{
  NextPartition();
  return *this;
}

G4FermiIntegerPartition::Iterator G4FermiIntegerPartition::Iterator::operator++(int)
{
  auto copy = *this;
  NextPartition();
  return copy;
}

bool G4FermiIntegerPartition::Iterator::operator==(
  const G4FermiIntegerPartition::Iterator& other) const
{
  return partition_ == other.partition_;
}

bool G4FermiIntegerPartition::Iterator::operator!=(
  const G4FermiIntegerPartition::Iterator& other) const
{
  return partition_ != other.partition_;
}

G4FermiIntegerPartition::Iterator::Iterator(uint32_t number, uint32_t termsCount, uint32_t base)
  : partition_(termsCount, 0)
{
  // No possible partitions
  if (number < base * termsCount || termsCount == 0 || number == 0) {
    return;
  }

  std::fill(partition_.begin(), partition_.end(), base);
  partition_[0] = number - base * (termsCount - 1);
}

void G4FermiIntegerPartition::Iterator::NextPartition()
{
  uint32_t accumulated = 0;
  for (auto partitionLast = std::next(partition_.begin()); partitionLast != partition_.end();
       ++partitionLast)
  {
    if (partition_.front() >= *partitionLast + 2) {
      --partition_.front();
      ++(*partitionLast);

      auto newValue = *partitionLast;
      std::fill(std::next(partition_.begin()), partitionLast, newValue);
      partition_.front() +=
        accumulated - newValue * (std::distance(partition_.begin(), partitionLast) - 1);
      return;
    }
    accumulated += *partitionLast;
  }

  // last partition
  partition_.clear();
}
