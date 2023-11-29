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
//

#ifndef G4OctreeFinder_hh
#define G4OctreeFinder_hh 1
#include "globals.hh"
#include <map>
#include "G4Octree.hh"
#include "G4Track.hh"
#include "G4ITType.hh"
#include "G4memory.hh"
#include "G4TrackList.hh"

#undef DEBUG
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4VFinder
{
public:
    G4VFinder() = default;
    virtual ~G4VFinder() = default;
    virtual void Clear() = 0;
    virtual void SetVerboseLevel(G4int level) = 0;
    virtual G4int GetVerboseLevel() = 0;
    virtual G4ITType GetITType() = 0;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef _Extractor_
#define _Extractor_
template<typename CONTAINER>
class Extractor
{
public:
    const G4ThreeVector operator () (typename CONTAINER::iterator it)
    {
	  return (*it)->GetPosition();
    }
    std::function<G4bool(const std::pair<typename CONTAINER::iterator,G4double>,
                         const std::pair<typename CONTAINER::iterator,G4double>)>
    compareInterval = [](const std::pair<typename CONTAINER::iterator,G4double>& iter1,
                         const std::pair<typename CONTAINER::iterator,G4double>& iter2)
    -> G4bool
    {
        return (std::get<1>(iter1) < std::get<1>(iter2));
    };
    
};
#endif
template<class T,typename CONTAINER>
class G4OctreeFinder: public G4VFinder
{
    using Octree = G4Octree<typename CONTAINER::iterator,
    Extractor<CONTAINER> >;
    using OctreeHandle = G4shared_ptr<Octree>;
    using TreeMap = std::map<int, OctreeHandle>;
    
private:
    static G4ThreadLocal G4OctreeFinder* fInstance;
    G4OctreeFinder();
    int fVerbose;
    G4bool fIsOctreeUsed;
    G4bool fIsOctreeBuit;
    Extractor<CONTAINER> fExtractor;
    TreeMap fTreeMap;
    OctreeHandle fTree;
public:
    static G4OctreeFinder * Instance();
    
    void SetOctreeUsed(G4bool used);
    G4bool IsOctreeUsed() const;
    
    void SetOctreeBuilt(G4bool used);
    G4bool IsOctreeBuilt() const;
    
    ~G4OctreeFinder() override;
    void Clear() override;

    void SetVerboseLevel(G4int level) override
    {
        fVerbose = level;
    }

    G4int GetVerboseLevel() override
    {
        return fVerbose;
    }

    G4ITType GetITType() override
    {
        return T::ITType();
    }
    void BuildTreeMap(const std::map<G4int,CONTAINER*>& listMap);
    void FindNearestInRange(const G4Track& track,
                            const int& key,
                            G4double R,
                            std::vector<std::pair<typename
                            CONTAINER::iterator,G4double>>& result,
                            G4bool isSort = false) const;

    void FindNearest(const G4Track& track,
                            const int& key,
                            G4double R,
                            std::vector<std::pair<typename
                            CONTAINER::iterator,G4double>>& result,
                            G4bool isSort = false) const;
    
    void FindNearestInRange(const G4ThreeVector& position,
                            const G4int& key,
                            G4double R,
                            std::vector<std::pair<typename
                            CONTAINER::iterator,G4double>>& result,
                            G4bool isSort = false) const;

    void FindNearestInRange(const G4ThreeVector& /*from this point*/,
                            G4double R,
                            std::vector<std::pair<
                            typename CONTAINER::iterator,G4double> >&
                            result,
                            G4bool isSorted) const;
};

#include "G4OctreeFinder.icc"
#endif