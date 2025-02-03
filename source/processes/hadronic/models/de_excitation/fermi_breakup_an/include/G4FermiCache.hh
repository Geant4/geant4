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

#ifndef G4FERMICACHE_HH
#define G4FERMICACHE_HH

#include "G4FermiDataTypes.hh"
#include "G4FermiLogger.hh"

#include <cassert>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>

template<typename Key, typename Value>
class G4FermiLFUCache : public G4FermiVCache<Key, Value>
{
  public:
    G4FermiLFUCache(size_t maxSize) : frequencyStorage_(), cache_(), maxSize_(maxSize)
    {
      FERMI_ASSERT_MSG(maxSize_ > 0, "Cahce size must be positive");
    }

    std::shared_ptr<Value> Insert(const Key& key, Value&& value) override
    {
      if (auto it = cache_.find(key); it != cache_.end()) {
        frequencyStorage_.erase(it->second.second);
        it->second.second = frequencyStorage_.insert({1, key});
        return it->second.first;
      }

      if (cache_.size() == maxSize_) {
        auto lowestIt = frequencyStorage_.begin();
        cache_.erase(lowestIt->second);
        frequencyStorage_.erase(lowestIt);
      }

      auto it = frequencyStorage_.insert({1, key});
      auto [insertedIt, _] =
        cache_.emplace(key, std::make_pair(std::make_shared<Value>(std::move(value)), it));
      return insertedIt->second.first;
    }

    std::shared_ptr<Value> Get(const Key& key) override
    {
      if (auto it = cache_.find(key); it != cache_.end()) {
        auto& freqIter = it->second.second;
        const auto newCount = freqIter->first + 1;
        frequencyStorage_.erase(freqIter);
        it->second.second = frequencyStorage_.insert({newCount, key});
        return it->second.first;
      }

      return nullptr;
    }

  private:
    using G4FermiFrequencyStorage = std::multimap<std::size_t, Key>;
    using G4FermiFrequencyIter = typename G4FermiFrequencyStorage::iterator;

    G4FermiFrequencyStorage frequencyStorage_;
    std::unordered_map<Key, std::pair<std::shared_ptr<Value>, G4FermiFrequencyIter>> cache_;
    std::size_t maxSize_;
};

template<typename Key, typename Value>
class G4FermiSimpleCache : public G4FermiVCache<Key, Value>
{
  public:
    G4FermiSimpleCache() = default;

    std::shared_ptr<Value> Insert(const Key& key, Value&& value) override
    {
      cachedItem_ = {key, std::make_shared<Value>(std::move(value))};
      return cachedItem_->second;
    }

    std::shared_ptr<Value> Get(const Key& key) override
    {
      if (cachedItem_.has_value() && cachedItem_->first == key) {
        return cachedItem_->second;
      }

      return nullptr;
    }

  private:
    std::optional<std::pair<Key, std::shared_ptr<Value>>> cachedItem_;
};

#endif  // G4FERMICACHE_HH
