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

#ifndef G4FERMIINTEGERPARTITION_HH
#define G4FERMIINTEGERPARTITION_HH

#include <cstddef>
#include <cstdint>
#include <vector>

namespace fbu
{

using G4FermiPartition = std::vector<uint32_t>;

class G4FermiIntegerPartition
{
  public:
    class Iterator;

    Iterator begin() const;

    Iterator end() const;

    G4FermiIntegerPartition(uint32_t number, uint32_t termsCount, uint32_t base = 1);

  private:
    uint32_t number_;
    uint32_t termsCount_;
    uint32_t base_;
};

class G4FermiIntegerPartition::Iterator
{
  public:
    friend class G4FermiIntegerPartition;

    using difference_type = int64_t;
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

    bool operator==(const Iterator& other) const;

    bool operator!=(const Iterator& other) const;

  private:
    // represents end partition
    Iterator() = default;

    Iterator(uint32_t number, uint32_t termsCount, uint32_t base);

    void NextPartition();

    G4FermiPartition partition_;
};

}  // namespace fbu

#endif  // G4FERMIINTEGERPARTITION_HH
