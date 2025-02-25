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

// Class template for accumulable maps handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 19/07/2024

#ifndef G4AccMap_h
#define G4AccMap_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

#include <map>

template <class Key,
          class T,
          class Compare = std::less<Key>,
          class Allocator = std::allocator<std::pair<const Key, T>>>
class G4AccMap : public G4VAccumulable
{
  public:
    // ctors to be supported
    // (https://en.cppreference.com/w/cpp/container/map/map)
    //
    // Default constructor (1) - Constructs an empty container.
    // map();
    //
    // Constructor (2) - Constructs an empty container.
    // explicit map(const Compare& comp,
    //              const Allocator& alloc = Allocator());
    //
    // Constructor (3) - Constructs an empty container.
    // explicit map(const Allocator& alloc);
    //
    // 4,5) Constructs the container with the contents of the range [first, last).
    // - skipped
    //
    // Constructor (6) - Copy constructor.
    // Constructs the container with the copy of the contents of other.
    // map( const map& other );
    //
    // Constructor (7) - Copy constructor with provided allocator
    // map( const map& other, const Allocator& alloc );
    //
    // Constructor (8) - Move constructor.
    // Constructs the container with the contents of other using move semantics.
    // map(map&& other);
    //
    // Constructor (9) - Move constructor with provided allocator
    // Constructs the container with the contents of other using move semantics.
    // map(map&& other, const Allocator& alloc);
    //
    // Constructor (10) - Initializer-list constructor.
    // Constructs the container with the contents of the initializer list init
    // map(std::initializer_list<value_type> init,
    //     const Compare& comp = Compare(),
    //     const Allocator& alloc = Allocator());
    //
    // Constructor (11) - Initializer-list constructor.
    // Constructs the container with the contents of the initializer list init
    // map(std::initializer_list<value_type> init,
    //     const Allocator& alloc);
    // - skipped
    //
    // (since C++20)
    // (12,13) Constructs the container with the contents of rg.
    // - skipped

    // Default constructor (1)
    // Constructs an empty container with all defaults.
    G4AccMap(const G4String& name = "",
             G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (2)
    // Constructs an empty container with the given compare.
    G4AccMap(const Compare& comp,
             G4MergeMode mergeMode = G4MergeMode::kAddition,
             const Allocator& alloc = Allocator());

    // Constructor (2) with name
    // Constructs an empty container with the given compare with name
    G4AccMap(const G4String& name,
             const Compare& comp,
             G4MergeMode mergeMode = G4MergeMode::kAddition,
             const Allocator& alloc = Allocator());

    // Constructor (3)
    // Constructs an empty container with the given allocator alloc.
    G4AccMap(const Allocator& alloc,
             G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (3) with name
    // Constructs an empty container with the given allocator alloc and name
    G4AccMap(const G4String& name,
             const Allocator& alloc,
             G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (10)
    // Constructs the container with the contents of the initializer list init.
    G4AccMap(std::initializer_list<std::pair<const Key,T>> init,
             G4MergeMode mergeMode = G4MergeMode::kAddition,
             const Compare& comp = Compare(),
             const Allocator& alloc = Allocator());

    // Constructor (10) with name
    // Constructs the container with the contents of the initializer list init.
    G4AccMap(const G4String& name,
             std::initializer_list<std::pair<const Key,T>> init,
             G4MergeMode mergeMode = G4MergeMode::kAddition,
             const Compare& comp = Compare(),
             const Allocator& alloc = Allocator());

    // Copy constructor
    G4AccMap(const G4AccMap& rhs) = default;
    G4AccMap(const G4AccMap& rhs, const Allocator& allocator);
    // Move constructor
    G4AccMap(G4AccMap&& rhs) = default;
    G4AccMap(G4AccMap&& rhs, const Allocator& allocator);

    // Destructor
    ~G4AccMap() override = default;

    // std::map functions, operators
    // operator []
    inline T& operator[](const Key& key) { return fMap[key]; }
    inline T& operator[](Key&& key) { return fMap[std::move(key)]; }
    // at
    inline T& at(const Key& key) { return fMap[key]; }
    inline const T& at(const Key& key ) const { return fMap[key]; }
    // size
    inline typename std::map<Key, T, Compare, Allocator>::size_type size() const { return fMap.size(); }
    // begin, cbegin
    inline typename std::map<Key, T, Compare, Allocator>::iterator begin() { return fMap.begin(); }
    inline typename std::map<Key, T, Compare, Allocator>::const_iterator begin() const { return fMap.begin(); }
    inline typename std::map<Key, T, Compare, Allocator>::const_iterator cbegin() const { return fMap.cbegin(); }
    // end, cend
    inline typename std::map<Key, T, Compare, Allocator>::iterator end() { return fMap.end(); }
    inline typename std::map<Key, T, Compare, Allocator>::const_iterator end() const { return fMap.end(); }
    inline typename std::map<Key, T, Compare, Allocator>::const_iterator cend() const { return fMap.cend(); }
    // clear
    inline void clear() { fMap.clear(); }
    // insert
    inline std::pair<typename std::map<Key, T, Compare, Allocator>::iterator, bool> insert(const T& value) { return fMap.insert(value); }
    template< class P >
    inline std::pair<typename std::map<Key, T, Compare, Allocator>::iterator, bool> insert( P&& value ) { return fMap.insert(std::move(value)); }
    inline std::pair<typename std::map<Key, T, Compare, Allocator>::iterator, bool> insert( T&& value ) { return fMap.insert(std::move(value)); }   
    // find
    inline typename std::map<Key, T, Compare, Allocator>::iterator find( const Key& key ) { return fMap.find(key); }
    inline typename std::map<Key, T, Compare, Allocator>::const_iterator find( const Key& key ) const { return fMap.find(key); }

    // Methods
    void Merge(const G4VAccumulable& other) final;
    void Reset() final;
    void Print(G4PrintOptions options = G4PrintOptions()) const final;
    void SetMergeMode(G4MergeMode value) final;
    void SetInitValue(const T& value);

    // Get methods
    G4AccType GetType() const final { return G4AccType::kMap; }
    std::map<Key, T, Compare, Allocator>& GetMap() { return fMap; }
    const std::map<Key, T, Compare, Allocator>& GetMap() const { return fMap; }

  private:
    // Data members
    std::map<Key, T, Compare, Allocator> fMap {};
    T fInitValue = 0;
    G4MergeFunction<T> fMergeFunction;
 };

// inline functions

#include "G4AccMap.icc"

#endif
