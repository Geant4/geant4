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

// The common implementation of analysis manager classes.

// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 07/09/2015

#ifndef G4AccumulableManager_h
#define G4AccumulableManager_h 1

#include "G4AccValue.hh"
#include "G4AccType.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4AccumulableManager;
class G4VAccumulable;
template <class T>
class G4ThreadLocalSingleton;
template <class T, std::size_t N>
class G4AccArray;
template <class Key, class T, class Compare, class Allocator>
class G4AccMap;
template <class Key, class T, class Hash, class KeyEqual, class Allocator>
class G4AccUnorderedMap;
template <class T, class Allocator>
class G4AccVector;

class G4AccumulableManager
{
  friend class G4ThreadLocalSingleton<G4AccumulableManager>;

  public:
    virtual ~G4AccumulableManager();

    // Static methods
    static G4AccumulableManager* Instance();

    // Methods

    // Create accumulables
    //
    template <typename T>
    G4AccValue<T>*
    CreateAccValue(const G4String& name, T value,
                   G4MergeMode mergeMode = G4MergeMode::kAddition);

    template <typename T>
    G4AccValue<T>*
    CreateAccValue(T value,
                   G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Deprecated function using the old accumulable value name
    //

    template <typename T>
    [[deprecated("Use `G4AccumulableManager::CreateAccValue<T>` instead")]]
    G4AccValue<T>*
    CreateAccumulable(const G4String& name, T value,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    template <typename T>
    [[deprecated("Use `G4AccumulableManager::CreateAccValue<T>` instead")]]
    G4AccValue<T>*
    CreateAccumulable(T value,
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Register existing accumulables
    //
    template <typename T>
    G4bool Register(G4AccValue<T>& accumulable);

    template <class T, std::size_t N>
    G4bool Register(G4AccArray<T, N>& accumulableArray);

    template <class Key, class T, class Compare, class Allocator>
    G4bool Register(G4AccMap<Key, T, Compare, Allocator>& accumulableMap);

    template <class Key, class T, class Hash, class KeyEqual, class Allocator>
    G4bool Register(G4AccUnorderedMap<Key, T, Hash, KeyEqual, Allocator>& accumulableUnorderedMap);

    template <class T, class Allocator>
    G4bool Register(G4AccVector<T, Allocator>& accumulableVector);

    // user defined accumulable
    G4bool Register(G4VAccumulable* accumulable);

    // Deprecated functions with long name
    //
    template <typename T>
    [[deprecated("Use `G4AccumulableManager::Register` instead")]]
    G4bool RegisterAccumulable(G4AccValue<T>& accumulable);
    // user defined accumulable
    [[deprecated("Use `G4AccumulableManager::Register` instead")]]
    G4bool RegisterAccumulable(G4VAccumulable* accumulable);

    // Access registered accumulables
    //
    // Via name
    // templated accumulable
    //
    template <typename T>
    G4AccValue<T>* GetAccValue(const G4String& name, G4bool warn = true) const;

    template <class T, std::size_t N>
    G4AccArray<T, N>*
    GetAccArray(const G4String& name, G4bool warn = true) const;

    template <class Key, class T, class Compare = std::less<Key>,
              class Allocator  = std::allocator<std::pair<const Key, T>>>
    G4AccMap<Key, T, Compare, Allocator>*
    GetAccMap(const G4String& name, G4bool warn = true) const;

    template <class Key, class T, class Hash = std::hash<Key>, class KeyEqual = std::equal_to<Key>,
              class Allocator  = std::allocator<std::pair<const Key, T>>>
    G4AccUnorderedMap<Key, T, Hash, KeyEqual, Allocator>*
    GetAccUnorderedMap(const G4String& name, G4bool warn = true) const;

    template <class T, class Allocator = std::allocator<T>>
    G4AccVector<T, Allocator>*
    GetAccVector(const G4String& name, G4bool warn = true) const;

    // user defined accumulable
    G4VAccumulable*  GetAccumulable(const G4String& name, G4bool warn = true) const;

    // Deprecated function using the old accumulable value name
    template <typename T>
    [[deprecated("Use `G4AccumulableManager::GetAccValue<T>` instead")]]
    G4AccValue<T>* GetAccumulable(const G4String& name, G4bool warn = true) const;

    // Via id (in the order of registering)
    // templated accumulable
    //
    template <typename T>
    G4AccValue<T>*  GetAccValue(G4int id, G4bool warn = true) const;

    template <class T, std::size_t N>
    G4AccArray<T, N>* GetAccArray(G4int id, G4bool warn = true) const;

    template <class Key, class T, class Compare = std::less<Key>,
              class Allocator = std::allocator<std::pair<const Key, T>>>
    G4AccMap<Key, T, Compare, Allocator>*
    GetAccMap(G4int id, G4bool warn = true) const;

    template <class Key, class T, class Hash = std::hash<Key>, class KeyEqual = std::equal_to<Key>,
              class Allocator = std::allocator<std::pair<const Key, T>>>
    G4AccUnorderedMap<Key, T, Hash, KeyEqual, Allocator>*
    GetAccUnorderedMap(G4int id, G4bool warn = true) const;

    template <class T, class Allocator = std::allocator<T>>
    G4AccVector<T, Allocator>*
    GetAccVector(G4int id, G4bool warn = true) const;

    // user defined accumulable
    G4VAccumulable*  GetAccumulable(G4int id, G4bool warn = true) const;

    // Deprecated function using the old accumulable value name
    template <typename T>
    [[deprecated("Use `G4AccumulableManager::GetAccValue<T>` instead")]]
    G4AccValue<T>*  GetAccumulable(G4int id, G4bool warn = true) const;

    G4int GetNofAccumulables() const;

    // Via vector iterators
    std::vector<G4VAccumulable*>::iterator Begin();
    std::vector<G4VAccumulable*>::iterator End();
    std::vector<G4VAccumulable*>::const_iterator BeginConst() const;
    std::vector<G4VAccumulable*>::const_iterator EndConst() const;

    // Methods applied to all accumulables
    void Merge();
    void Reset();
    void Print(G4PrintOptions options = G4PrintOptions()) const;

    // Print a selected range
    void Print(G4int startId, G4int count,
               G4PrintOptions options = G4PrintOptions()) const;
    void Print(std::vector<G4VAccumulable*>::iterator startIt, std::size_t count,
               G4PrintOptions options = G4PrintOptions()) const;
    void Print(std::vector<G4VAccumulable*>::iterator startIt,
               std::vector<G4VAccumulable*>::iterator endIt,
               G4PrintOptions options = G4PrintOptions()) const;

    // Verbose level
    void SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;

  private:
    // Hide singleton ctor
    G4AccumulableManager();

    // Methods
    // Generate generic accumulable name: accumulableN, where N is the actual number of accumulables
    G4String GenerateName() const;
    // Check if a name is already used in a map and print a warning
    G4bool CheckName(const G4String& name, const G4String& where) const;
    G4bool CheckType(G4VAccumulable* accumulable, G4AccType type, G4bool warn) const;

    template <typename T>
    G4AccValue<T>*  GetAccumulable(G4VAccumulable* accumulable, G4bool warn) const;

    template <typename T, std::size_t N>
    G4AccArray<T, N>* GetAccArray(G4VAccumulable* accumulable, G4bool warn) const;

    template <typename Key, typename T, typename Compare = std::less<Key>,
              typename Allocator = std::allocator<std::pair<const Key, T>>>
    G4AccMap<Key, T, Compare, Allocator>*
    GetAccMap(G4VAccumulable* accumulable, G4bool warn = true) const;

    template <typename Key, typename T, typename Hash = std::hash<Key>, typename KeyEqual = std::equal_to<Key>,
              typename Allocator = std::allocator<std::pair<const Key, T>>>
    G4AccUnorderedMap<Key, T, Hash, KeyEqual, Allocator>*
    GetAccUnorderedMap(G4VAccumulable* accumulable, G4bool warn = true) const;

    template <typename T, typename Allocator = std::allocator<T>>
    G4AccVector<T, Allocator>*
    GetAccVector(G4VAccumulable* accumulable, G4bool warn) const;

    // Constants
    const G4String kBaseName = "accumulable";

    // Static data members
    inline static G4AccumulableManager* fgMasterInstance { nullptr };

    // Data members
    std::vector<G4VAccumulable*>        fVector;
    std::map<G4String, G4VAccumulable*> fMap;
    std::vector<G4VAccumulable*>        fAccumulablesToDelete;
 };

#include "G4AccumulableManager.icc"

#endif

