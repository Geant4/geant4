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
//

#ifndef G4HumanPhantomHit_h
#define G4HumanPhantomHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh" // FOR MT

class G4HumanPhantomHit : public G4VHit
{
public:

  G4HumanPhantomHit();
  ~G4HumanPhantomHit();
  G4HumanPhantomHit(const G4HumanPhantomHit&);
  const G4HumanPhantomHit& operator=(const G4HumanPhantomHit&);
  G4bool operator==(const G4HumanPhantomHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

public:  
  void SetBodyPartID  (G4String bodyPartName) { bodyPartID = bodyPartName;};
  void SetEdep (G4double de) { edep = de; };
      
  G4String GetBodyPartID() { return bodyPartID; };
  G4double GetEdep()    { return edep; };      
      
private:
  G4String bodyPartID;
  G4double      edep;    
};

typedef G4THitsCollection<G4HumanPhantomHit> G4HumanPhantomHitsCollection;

extern G4ThreadLocal G4Allocator<G4HumanPhantomHit>* G4HumanPhantomHitAllocator;

inline void* G4HumanPhantomHit::operator new(size_t)
{
  if(!G4HumanPhantomHitAllocator)
      G4HumanPhantomHitAllocator = new G4Allocator<G4HumanPhantomHit>;
  return (void *) G4HumanPhantomHitAllocator->MallocSingle();
}

inline void G4HumanPhantomHit::operator delete(void *aHit)
{
  G4HumanPhantomHitAllocator -> FreeSingle((G4HumanPhantomHit*) aHit);
}

#endif
