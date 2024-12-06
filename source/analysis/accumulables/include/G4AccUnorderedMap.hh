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

// Class template for accumulable unordered_maps handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 19/07/2024

#ifndef G4AccUnorderedMap_h
#define G4AccUnorderedMap_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

#include <unordered_map>

template <class Key,
          class T,
          class Hash = std::hash<Key>,
          class KeyEqual = std::equal_to<Key>,
          class Allocator = std::allocator<std::pair<const Key, T>>>
class G4AccUnorderedMap : public G4VAccumulable
{
  public:
    // ctors to be supported
    // (https://en.cppreference.com/w/cpp/container/unordered_map/unordered_map)
    //
    // Default constructor (1) - Constructs an empty container.
    // unordered_map();
    //
    // Constructor (2) - Constructs an empty container.
    // explicit unordered_map( size_type bucket_count,
    //                         const Hash& hash = Hash(),
    //                         const KeyEqual& equal = KeyEqual(),
    //                         const Allocator& alloc = Allocator() );
    //
    // Constructor (3) - Constructs an empty container.
    // unordered_map( size_type bucket_count,
    //                const Allocator& alloc )
    //   : unordered_map(bucket_count, Hash(), KeyEqual(), alloc) {}
    //
    // Constructor (4)  (since C++14)
    // Constructs empty container.
    // unordered_map( size_type bucket_count,
    //                const Hash& hash,
    //                const Allocator& alloc )
    //     : unordered_map(bucket_count, hash, KeyEqual(), alloc) {}
    //
    // Constructor  (5) (since C++11)
    // Constructs empty container.
    // explicit unordered_map( const Allocator& alloc );
    //
    // (6,7,8) (since C++11)
    // Constructs the container with the contents of the range
    // - skipped
    //
    // (9) (since C++11)  - Copy constructor.
    // unordered_map( const unordered_map& other );
    //
    // (10) (since C++11) - Copy constructor with provided allocator
    // unordered_map( const unordered_map& other, const Allocator& alloc );
    //
    // (11)  (since C++11) - Move constructor.
    // unordered_map( unordered_map&& other );
    //
    // (12)  (since C++11) - Move constructor with provided allocator
    // unordered_map( unordered_map&& other, const Allocator& alloc );
    //
    // (13)  (since C++11) - Initializer-list constructor.
    // Constructs the container with the contents of the initializer list init
    // unordered_map( std::initializer_list<value_type> init,
    //                size_type bucket_count = /* implementation-defined */,
    //                const Hash& hash = Hash(),
    //                const KeyEqual& equal = KeyEqual(),
    //                const Allocator& alloc = Allocator() );
    //
    // (14)  (since C++11) - Initializer-list constructor.
    // Constructs the container with the contents of the initializer list init
    // unordered_map( std::initializer_list<value_type> init,
    //                size_type bucket_count,
    //                const Allocator& alloc )
    //     : unordered_map(init, bucket_count,
    //                     Hash(), KeyEqual(), alloc) {}
    //
    // (15)  (since C++14) - Initializer-list constructor.
    // unordered_map( std::initializer_list<value_type> init,
    //                size_type bucket_count,
    //                const Hash& hash,
    //                const Allocator& alloc )
    //     : unordered_map(init, bucket_count,
    //                     hash, KeyEqual(), alloc) {}
    //
    // (16,17)  (since C++23)
    // Constructs the container with the contents of rg.
    // - skipped


    // Default constructor (1)
    // Constructs an empty container with all defaults.
    G4AccUnorderedMap(const G4String& name = "",
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (2)
    // Constructs an empty container with the given bucket_count
    G4AccUnorderedMap(std::size_t bucket_count,
                      G4MergeMode mergeMode = G4MergeMode::kAddition,
                      const Allocator& alloc = Allocator());

    // Constructor (2) with name
    // Constructs an empty container with the given bucket_count and name
    G4AccUnorderedMap(const G4String& name,
                      std::size_t bucket_count,
                      G4MergeMode mergeMode = G4MergeMode::kAddition,
                      const Allocator& alloc = Allocator());

    // Constructor (3)
    // Constructs an empty container with the given bucket_count and allocator
    G4AccUnorderedMap(std::size_t bucket_count,
                      const Allocator& alloc,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (3) with name
    // Constructs an empty container with the given bucket_count, allocator and name
    G4AccUnorderedMap(const G4String& name,
                      std::size_t bucket_count,
                      const Allocator& alloc,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);
    // Constructor (4)
    // Constructs an empty container with the given bucket_count, allocator and hash
    G4AccUnorderedMap(std::size_t bucket_count,
                      const Hash& hash,
                      const Allocator& alloc,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (4) with name
    // Constructs an empty container with the given bucket_count, allocator, hash and name
    G4AccUnorderedMap(const G4String& name,
                      std::size_t bucket_count,
                      const Hash& hash,
                      const Allocator& alloc,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);
    // Constructor (5)
    // Constructs an empty container with the given allocator alloc.
    G4AccUnorderedMap(const Allocator& alloc,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (5) with name
    // Constructs an empty container with the given allocator alloc and name
    G4AccUnorderedMap(const G4String& name,
                      const Allocator& alloc,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (13)
    // Constructs the container with the contents of the initializer list init.
    G4AccUnorderedMap(std::initializer_list<std::pair<const Key,T>> init,
                      G4MergeMode mergeMode = G4MergeMode::kAddition,
                      std::size_t bucket_count = 0,
                      const Hash& hash = Hash(),
                      const KeyEqual& equal = KeyEqual(),
                      const Allocator& alloc = Allocator() );

    // Constructor (13) with name
    // Constructs the container with the contents of the initializer list init.
    G4AccUnorderedMap(const G4String& name,
                      std::initializer_list<std::pair<const Key,T>> init,
                      G4MergeMode mergeMode = G4MergeMode::kAddition,
                      std::size_t bucket_count = 0,
                      const Hash& hash = Hash(),
                      const KeyEqual& equal = KeyEqual(),
                      const Allocator& alloc = Allocator() );

    // Copy constructor
    G4AccUnorderedMap(const G4AccUnorderedMap& rhs) = default;
    G4AccUnorderedMap(const G4AccUnorderedMap& rhs, const Allocator& allocator);
    // Move constructor
    G4AccUnorderedMap(G4AccUnorderedMap&& rhs) = default;
    G4AccUnorderedMap(G4AccUnorderedMap&& rhs, const Allocator& allocator);

    // Destructor
    ~G4AccUnorderedMap() override = default;

    // std::unordered_map functions, operators
    // operator []
    inline T& operator[](const Key& key) { return fUMap[key]; }
    inline T& operator[](Key&& key) { return fUMap[std::move(key)]; }
    // at
    inline T& at(const Key& key) { return fUMap[key]; }
    inline const T& at(const Key& key ) const { return fUMap[key]; }
    // size
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::size_type size() const { return fUMap.size(); }
    // begin, cbegin
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::iterator begin() { return fUMap.begin(); }
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::const_iterator begin() const { return fUMap.begin(); }
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::const_iterator cbegin() const { return fUMap.cbegin(); }
    // end, cend
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::iterator end() { return fUMap.end(); }
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::const_iterator end() const { return fUMap.end(); }
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::const_iterator cend() const { return fUMap.cend(); }
    // clear
    inline void clear() { fUMap.clear(); }
    // insert
    inline std::pair<typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::iterator, bool> insert(const T& value) { return fUMap.insert(value); }
    template< class P >
    inline std::pair<typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::iterator, bool> insert( P&& value ) { return fUMap.insert(std::move(value)); }
    inline std::pair<typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::iterator, bool> insert( T&& value ) { return fUMap.insert(std::move(value)); }   
    // find
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::iterator find( const Key& key ) { return fUMap.find(key); }
    inline typename std::unordered_map<Key, T, Hash, KeyEqual, Allocator>::const_iterator find( const Key& key ) const { return fUMap.find(key); }

    // Methods
    void Merge(const G4VAccumulable& other) final;
    void Reset() final;
    void Print(G4PrintOptions options = G4PrintOptions()) const final;
    void SetMergeMode(G4MergeMode value) final;
    void SetInitValue(const T& value);

    // Get methods
    G4AccType GetType() const final { return G4AccType::kUnorderedMap; }
    std::unordered_map<Key, T, Hash, KeyEqual, Allocator>& GetUnorderedMap() { return fUMap; }
    const std::unordered_map<Key, T, Hash, KeyEqual, Allocator>& GetUnorderedMap() const { return fUMap; }

  private:
    // Data members
    std::unordered_map<Key, T, Hash, KeyEqual, Allocator> fUMap;
    T fInitValue = 0;
    G4MergeFunction<T> fMergeFunction;
 };

// inline functions

#include "G4AccUnorderedMap.icc"

#endif
