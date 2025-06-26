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
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#include "G4FermiIntegerPartition.hh"

G4integerPartition::G4integerPartition(std::uint32_t number, std::uint32_t termsCount,
                                       std::uint32_t base)
  : number_(number), termsCount_(termsCount), base_(base)
{}

G4integerPartition::Iterator G4integerPartition::begin() const
{
  return {number_, termsCount_, base_};
}

G4integerPartition::Iterator G4integerPartition::end() const
{
  return {};
}

/////////////////////////////////// ITERATOR //////////////////////////////

G4integerPartition::Iterator::pointer G4integerPartition::Iterator::operator->() const
{
  return &partition_;
}

G4integerPartition::Iterator::reference G4integerPartition::Iterator::operator*() const
{
  return partition_;
}

G4integerPartition::Iterator& G4integerPartition::Iterator::operator++()
{
  NextPartition();
  return *this;
}

G4integerPartition::Iterator G4integerPartition::Iterator::operator++(int)
{
  auto copy = *this;
  NextPartition();
  return copy;
}

G4bool G4integerPartition::Iterator::operator==(const G4integerPartition::Iterator& other) const
{
  return partition_ == other.partition_;
}

G4bool G4integerPartition::Iterator::operator!=(const G4integerPartition::Iterator& other) const
{
  return partition_ != other.partition_;
}

G4integerPartition::Iterator::Iterator(std::uint32_t number, std::uint32_t termsCount,
                                       std::uint32_t base)
  : partition_(termsCount, 0)
{
  // No possible partitions
  if (number < base * termsCount || termsCount == 0 || number == 0) {
    return;
  }

  std::fill(partition_.begin(), partition_.end(), base);
  partition_[0] = number - base * (termsCount - 1);
}

void G4integerPartition::Iterator::NextPartition()
{
  std::uint32_t accumulated = 0;
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
