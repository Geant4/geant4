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
/// \file exoticphysics/phonon/include/XAlminumElectrodeHit.hh
/// \brief Definition of the XAlminumElectrodeHit class
//
// $Id$
//
#ifndef XAlminumElectrodeHit_h
#define XAlminumElectrodeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class XAlminumElectrodeHit : public G4VHit
{
  public:

      XAlminumElectrodeHit();
      virtual ~XAlminumElectrodeHit();
      XAlminumElectrodeHit(const XAlminumElectrodeHit &right);
      const XAlminumElectrodeHit& operator=(const XAlminumElectrodeHit &right);
      int operator==(const XAlminumElectrodeHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
      G4double time;
      G4double fEdep;
      G4ThreeVector fLocalPos;
      G4ThreeVector fWorldPos;

  public:
      inline void SetTime(G4double t) { time = t; }
      inline G4double GetTime() const { return time; }
      inline void SetEDep(G4double e) { fEdep = e; }
      inline G4double GetEDep() const { return fEdep; }
      inline void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
      inline G4ThreeVector GetLocalPos() const { return fLocalPos; }
      inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
      inline G4ThreeVector GetWorldPos() const { return fWorldPos; }
};

typedef G4THitsCollection<XAlminumElectrodeHit> XAlminumElectrodeHitsCollection;

extern G4Allocator<XAlminumElectrodeHit> XAlminumElectrodeHitAllocator;

inline void* XAlminumElectrodeHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)XAlminumElectrodeHitAllocator.MallocSingle();
  return aHit;
}

inline void XAlminumElectrodeHit::operator delete(void* aHit)
{
  XAlminumElectrodeHitAllocator.FreeSingle((XAlminumElectrodeHit*) aHit);
}

#endif


