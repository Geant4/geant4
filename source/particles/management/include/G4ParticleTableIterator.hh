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
// G4ParticleTableIterator

// Authors: G.Cosmo, 2 December 1995 - Design, based on object model
//          H.Kurashige, 28 October 1999 - First implementation
// --------------------------------------------------------------------
#ifndef G4ParticleTableIterator_hh
#define G4ParticleTableIterator_hh 1

#include "G4ParticleDefinition.hh"

#include <map>

template<class K, class V>
class G4ParticleTableIterator
{
  public:
    using Map = std::map<K, V, std::less<K>>;

    G4ParticleTableIterator(Map& adict) : it(adict.begin()), mydict(&adict) {}

    G4bool operator++()
    {
      if (!defined) return false;
      ++it;
      return static_cast<bool>(it != mydict->end());
    }

    G4bool operator()()
    {
      if (!defined) {
        defined = true;
        it = mydict->begin();
      }
      else {
        ++it;
      }
      if (it == mydict->end()) return false;
      if (skipIons) {
        while ((static_cast<G4ParticleDefinition*>((*it).second))->IsGeneralIon())
        {  // Loop checking, 09.08.2015, K.Kurashige
          ++it;
          if (it == mydict->end()) return false;
        }
      }
      return true;
    }

    void reset(G4bool ifSkipIon = true)
    {
      defined = false;
      skipIons = ifSkipIon;
    }

    K* key() const { return &((*it).first); }

    V value() const { return (*it).second; }

  private:
    typename Map::iterator it;
    Map* mydict = nullptr;
    G4bool defined = false;
    G4bool skipIons = true;
};

#endif
