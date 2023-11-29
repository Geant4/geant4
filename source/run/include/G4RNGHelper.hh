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
// G4TemplateRNGHelper
//
// Class description:
//
// Helper class for RNG Engine seeds.
// Used in MT builds to guarantee event reproducibility.
// The function of this class is to return a RNG Engine seed given its index.
// It is a simple templated container that allows to add seeds (AddOneSeed)
// and retrieve a seed (GetSeed) by index.
//
// The normal use is with G4RNGHelper where each element of the container
// represents a seed. To enforce strong-reproducibility the variant with
// the RNG Engine status file names is available.

// Author: A.Dotti (SLAC), 5 July 2013
// --------------------------------------------------------------------
#ifndef G4TemplateRNGHelper_hh
#define G4TemplateRNGHelper_hh 1

#include "globals.hh"

#include <queue>
#include <vector>

template<class T>
class G4TemplateRNGHelper
{
  public:
    using SeedsQueue = std::vector<T>;
    using SeedsQueueSize_type = typename SeedsQueue::size_type;

  public:
    // The container is modeled as a (shared) singleton
    static G4TemplateRNGHelper<T>* GetInstance();
    static G4TemplateRNGHelper<T>* GetInstanceIfExist();
    virtual ~G4TemplateRNGHelper();

    // Returns seed given id
    virtual const T GetSeed(const G4int& sdId)
    {
      G4int seedId = sdId - 2 * offset;
      if (seedId < static_cast<G4int>(seeds.size())) {
        T& seed = seeds[seedId];
        return seed;
      }
      G4ExceptionDescription msg;
      msg << "No seed number " << seedId << "(" << seeds.size() << " available)\n"
          << " Original seed number " << sdId << " filled so far " << offset;
      G4Exception("G4RNGHelper::GetSeed", "Run0115", FatalException, msg);
      return T();
    }

    // Adds one seed to the collection
    void AddOneSeed(const T& seed) { seeds.push_back(seed); }

    // Fills N primary seed pairs
    void Fill(G4double* dbl, G4int nev, G4int nev_tot, G4int nrpe)
    {
      seeds.clear();
      for (G4int i = 0; i < nrpe * nev; ++i) {
        seeds.push_back((G4long)(100000000L * dbl[i]));
      }
      offset = 0;
      nev_filled = nev;
      nev_total = nev_tot;
      nRandParEvent = nrpe;
    }

    void Refill(G4double* dbl, G4int nev)
    {
      if (nev == 0) return;
      seeds.clear();
      for (G4int i = 0; i < nRandParEvent * nev; ++i) {
        seeds.push_back((G4long)(100000000L * dbl[i]));
      }
      offset += nev_filled;
      nev_filled = nev;
    }

    // Number of available seeds
    const SeedsQueueSize_type GetNumberSeeds() const { return seeds.size(); }

    // Empty the seeds container
    virtual void Clear() { seeds.clear(); }

  protected:
    SeedsQueue seeds;
    // Note: following numbers are number of events.
    //       seeds are generated for nRandParEvent times n_event
    G4int offset = 0;
    G4int nev_filled = 0;
    G4int nev_total = 0;
    G4int nRandParEvent = 0;

  private:
    G4TemplateRNGHelper() = default;

  private:
    static G4TemplateRNGHelper<T>* instance;
};

using G4RNGHelper = G4TemplateRNGHelper<G4long>;
using G4StringRNGHelper = G4TemplateRNGHelper<G4String>;
using G4SeedsQueue = std::queue<G4long>;

#endif
