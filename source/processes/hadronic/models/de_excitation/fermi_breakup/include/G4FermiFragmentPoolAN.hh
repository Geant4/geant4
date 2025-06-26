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

#ifndef G4FERMIFRAGMENTPOOLAN_HH
#define G4FERMIFRAGMENTPOOLAN_HH

#include "G4FermiDataTypes.hh"
#include "G4VFermiFragmentAN.hh"
#include "globals.hh"

class G4FermiFragmentPoolAN
{
  private:
    using Container = std::vector<const G4VFermiFragmentAN*>;

  public:
    class DefaultPoolANSource : private std::vector<G4VFermiFragmentAN*>
    {
      private:
        using PoolANContainer = std::vector<G4VFermiFragmentAN*>;

      public:
        DefaultPoolANSource();

        void Initialize();

        using PoolANContainer::begin;
        using PoolANContainer::cbegin;
        using PoolANContainer::cend;
        using PoolANContainer::end;
    };

    class IteratorRange
    {
      public:
        using const_iterator = Container::const_iterator;

        IteratorRange(const_iterator begin, const_iterator end) : begin_(begin), end_(end) {}

        const_iterator begin() const { return begin_; }
        const_iterator end() const { return end_; }

      private:
        const_iterator begin_;
        const_iterator end_;
    };

    std::size_t Count(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber) const;

    std::size_t Count(G4FermiNucleiData nuclei) const
    {
      return Count(nuclei.atomicMass, nuclei.chargeNumber);
    }

    IteratorRange GetFragments(G4FermiAtomicMass atomicMass,
                               G4FermiChargeNumber chargeNumber) const;

    IteratorRange GetFragments(G4FermiNucleiData nuclei) const
    {
      return GetFragments(nuclei.atomicMass, nuclei.chargeNumber);
    }

    template<typename DataSource>
    void Initialize(const DataSource& dataSource)
    {
      Initialize(dataSource.begin(), dataSource.end());
    }

    template<typename Iter>
    void Initialize(Iter begin, Iter end)
    {
      fragments_.clear();
      static_assert(
        std::is_same_v<std::remove_const_t<typename Iter::value_type>, G4VFermiFragmentAN*>,
        "invalid iterator");
      for (auto it = begin; it != end; ++it) {
        AddFragment(**it);
      }
    }

    void AddFragment(const G4VFermiFragmentAN& fragment);

    static G4FermiFragmentPoolAN& Instance()
    {
      static G4FermiFragmentPoolAN pool;
      return pool;
    }

  private:
    G4FermiFragmentPoolAN();
  
    static inline const Container EmptyContainer_ = {};

    std::vector<Container> fragments_;
};

#endif  // G4FERMIFRAGMENTPOOL_HH
