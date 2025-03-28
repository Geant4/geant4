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

#ifndef G4intEGERPARTITION_HH
#define G4intEGERPARTITION_HH

#include <globals.hh>

using G4FermiPartition = std::vector<std::uint32_t>;

class G4integerPartition
{
  public:
    class Iterator;

    Iterator begin() const;

    Iterator end() const;

    G4integerPartition(std::uint32_t number, std::uint32_t termsCount, std::uint32_t base = 1);

  private:
    std::uint32_t number_;
    std::uint32_t termsCount_;
    std::uint32_t base_;
};

class G4integerPartition::Iterator
{
  public:
    friend class G4integerPartition;

    using difference_type = std::int64_t;
    using value_type = G4FermiPartition;
    using reference = const G4FermiPartition&;
    using pointer = const G4FermiPartition*;
    using iterator_category = std::forward_iterator_tag;

    Iterator(const Iterator&) = default;

    Iterator& operator=(const Iterator&) = default;

    pointer operator->() const;

    reference operator*() const;

    Iterator& operator++();

    Iterator operator++(int);

    G4bool operator==(const Iterator& other) const;

    G4bool operator!=(const Iterator& other) const;

  private:
    // represents end partition
    Iterator() = default;

    Iterator(std::uint32_t number, std::uint32_t termsCount, std::uint32_t base);

    void NextPartition();

    G4FermiPartition partition_;
};

#endif  // G4intEGERPARTITION_HH
