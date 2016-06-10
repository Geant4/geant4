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
/// \file exoticphysics/phonon/include/XAluminumElectrodeHit.hh
/// \brief Definition of the XAluminumElectrodeHit class
//
// $Id: XAluminumElectrodeHit.hh 84197 2014-10-10 14:33:03Z gcosmo $
//
// 20141008  Allocators must be thread-local, and must be pointers

#ifndef XAluminumElectrodeHit_h
#define XAluminumElectrodeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class XAluminumElectrodeHit : public G4VHit
{
  public:

      XAluminumElectrodeHit();
      virtual ~XAluminumElectrodeHit();
      XAluminumElectrodeHit(const XAluminumElectrodeHit &right);
      const XAluminumElectrodeHit& operator=(const XAluminumElectrodeHit &right);
      int operator==(const XAluminumElectrodeHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
      G4double fTime;
      G4double fEdep;
      G4ThreeVector fLocalPos;
      G4ThreeVector fWorldPos;

  public:
      inline void SetTime(G4double t) { fTime = t; }
      inline G4double GetTime() const { return fTime; }
      inline void SetEDep(G4double e) { fEdep = e; }
      inline G4double GetEDep() const { return fEdep; }
      inline void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
      inline G4ThreeVector GetLocalPos() const { return fLocalPos; }
      inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
      inline G4ThreeVector GetWorldPos() const { return fWorldPos; }
};

typedef G4THitsCollection<XAluminumElectrodeHit> XAluminumElectrodeHitsCollection;

extern G4ThreadLocal G4Allocator<XAluminumElectrodeHit>* XAluminumElectrodeHitAllocator;

inline void* XAluminumElectrodeHit::operator new(size_t)
{
  if (!XAluminumElectrodeHitAllocator)                        // Singleton
    XAluminumElectrodeHitAllocator = new G4Allocator<XAluminumElectrodeHit>;

  return (void*)XAluminumElectrodeHitAllocator->MallocSingle();
}

inline void XAluminumElectrodeHit::operator delete(void* aHit)
{
  XAluminumElectrodeHitAllocator->FreeSingle((XAluminumElectrodeHit*) aHit);
}

#endif


