// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02TrackerHitsCollection.ddl,v 1.3 1999/11/29 18:23:33 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef PersEx02TrackerHitsCollection_h
#define PersEx02TrackerHitsCollection_h 1

#include "HepODBMS/odbms/HepODBMS.h"
#include "HepODBMS/clustering/HepClusteringHint.h"

#include "G4PersistentTypes.hh"
#include "G4PVHit.hh"
#include "G4PVHitsCollection.hh"
#include "G4ThreeVector.hh"

#include "PersEx02TrackerHit.hh"

typedef d_Varray< d_Ref<PersEx02TrackerHit> > PersEx02TrackerHitsVArray;

class PersEx02TrackerHitsCollection
 : public G4PVHitsCollection
{
  public:
      PersEx02TrackerHitsCollection(G4String dName, G4String aName);
      ~PersEx02TrackerHitsCollection();

      void DrawAllHits();
      void PrintAllHits();

  public:
      static HepClusteringHint clustering;

  private:
      PersEx02TrackerHitsVArray hitsCollection;
      size_t elems;

  public:
      inline void Insert(d_Ref<PersEx02TrackerHit> pHit)
        {
          size_t curSize = hitsCollection.size();
          if (elems >= curSize)
            hitsCollection.resize(curSize + 1024);
          hitsCollection[elems++] = pHit;
        };
      inline void ShrinkToFit()
        { hitsCollection.resize(elems); };
      inline size_t size() const
      { return hitsCollection.size(); }
      inline HepRef(PersEx02TrackerHit) element(const size_t i) const
      { return hitsCollection[i]; }
      inline HepRef(PersEx02TrackerHit) operator[](const size_t i) const
      { return hitsCollection[i]; }

};

#endif

