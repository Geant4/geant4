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
// --------------------------------------------------------------
//

#ifndef SensitiveDetectorHit_h
#define SensitiveDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class SensitiveDetectorHit : public G4VHit
{
public:
    SensitiveDetectorHit();
    SensitiveDetectorHit(G4int z);
    virtual ~SensitiveDetectorHit();
    SensitiveDetectorHit(
                        const SensitiveDetectorHit &right);
    const SensitiveDetectorHit& operator=(
                        const SensitiveDetectorHit &right);
    int operator==(const SensitiveDetectorHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    
private:
    G4int fLayerID;
    G4ThreeVector fWorldPos;
    
public:
    inline void SetLayerID(G4int z) { fLayerID = z; }
    inline G4int GetLayerID() const { return fLayerID; }
    inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    inline G4ThreeVector GetWorldPos() const { return fWorldPos; }
};

typedef G4THitsCollection<SensitiveDetectorHit>
    SensitiveDetectorHitsCollection;

#ifdef G4MULTITHREADED
extern G4ThreadLocal G4Allocator<SensitiveDetectorHit>*
    SensitiveDetectorHitAllocator;
#else
extern G4Allocator<SensitiveDetectorHit>
    SensitiveDetectorHitAllocator;
#endif

inline void* SensitiveDetectorHit::operator new(size_t)
{
#ifdef G4MULTITHREADED
    if(!SensitiveDetectorHitAllocator)
        SensitiveDetectorHitAllocator =
        new G4Allocator<SensitiveDetectorHit>;
    return (void *) SensitiveDetectorHitAllocator->MallocSingle();
#else
    void* aHit;
    aHit = (void*)SensitiveDetectorHitAllocator.MallocSingle();
    return aHit;
#endif
}

inline void SensitiveDetectorHit::operator delete(void* aHit)
{
#ifdef G4MULTITHREADED
    SensitiveDetectorHitAllocator->FreeSingle((
                    SensitiveDetectorHit*) aHit);
#else
   SensitiveDetectorHitAllocator.FreeSingle((
                    SensitiveDetectorHit*) aHit);
#endif
}

#endif


