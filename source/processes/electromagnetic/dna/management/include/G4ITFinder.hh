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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4ITFinder_hh
#define G4ITFinder_hh 1

#include <G4AllITFinder.hh>

#include "globals.hh"
#include <map>
#include "G4KDTree.hh"
#include "G4KDTreeResult.hh"
#include "G4Track.hh"
#include "G4ITTrackHolder.hh"

class G4VITFinder
{
 public:
  G4VITFinder() = default;
  virtual ~G4VITFinder() = default;
  virtual void Clear()                      = 0;
  virtual void Push(G4Track* track)         = 0;
  virtual G4ITType GetITType()              = 0;
  virtual void UpdatePositionMap()          = 0;
 protected:
  G4int fVerbose = 0;
};

/**
 * Localize the nearest neighbor
 * For now, G4KDTree is used
 */

template <class T>  //, class SearcherT = G4KDTree>
class G4ITFinder : public G4VITFinder
{
 public:
  static G4ITFinder* Instance();
  ~G4ITFinder() override;
  void Clear() override;
  void SetVerboseLevel(G4int level) { fVerbose = level; }
  G4int GetVerboseLevel() { return fVerbose; }
  void Push(G4Track* track) override;
  G4ITType GetITType() override { return T::ITType(); }
  void UpdatePositionMap() override;
  G4KDTreeResultHandle FindNearestInRange(const T* point /*from this point*/,
                                          G4int key /*for this type*/,
                                          G4double /*range*/);
  G4KDTreeResultHandle FindNearest(const G4ThreeVector&,
                                   G4int key /*for this type*/);
  G4KDTreeResultHandle FindNearest(const T* /*from this point*/,
                                   G4int key /*for this type*/);
  G4KDTreeResultHandle FindNearestInRange(
    const G4ThreeVector& /*from this point*/, G4int key /*for this type*/,
    G4double /*range*/);
  G4KDTreeResultHandle FindNearest(const T* /*from this point*/,
                                   const T* /*for this type*/);
 private:
  static G4ThreadLocal G4ITFinder* fInstance;
  G4ITFinder();
  using TreeMap = std::map<G4int, G4KDTree*> ;
  TreeMap fTree;
};

#ifdef TEMPLATE
#  undef TEMPLATE
#endif

#define TEMPLATE template <class T>
#define G4ITMANAGER G4ITFinder<T>

#include "G4ITFinder.icc"

#undef TEMPLATE
#undef G4ITMANAGER

#endif
