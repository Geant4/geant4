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
// $Id: MyCalorimeterHit.hh,v 1.5 2006-06-29 21:46:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyCalorimeterHit_h
#define MyCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class MyCalorimeterHit : public G4VHit
{
  public:

      MyCalorimeterHit();
      ~MyCalorimeterHit();
      MyCalorimeterHit(const MyCalorimeterHit &right);
      MyCalorimeterHit(G4VPhysicalVolume* physVol);
      const MyCalorimeterHit& operator=(const MyCalorimeterHit &right);
      int operator==(const MyCalorimeterHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      G4VPhysicalVolume* pPhys;

  public:
      inline void SetEdep(G4double de)
      { edep = de; };
      inline void AddEdep(G4double de)
      { edep += de; };
      inline G4double GetEdep()
      { return edep; };
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; };
      inline G4ThreeVector GetPos()
      { return pos; };
      inline void SetRot(G4RotationMatrix rmat)
      { rot = rmat; };
      inline G4RotationMatrix GetRot()
      { return rot; };
      inline G4VPhysicalVolume* GetPPhys()
      { return pPhys; };

};

typedef G4THitsCollection<MyCalorimeterHit> MyCalorimeterHitsCollection;

extern G4Allocator<MyCalorimeterHit> MyCalorimeterHitAllocator;

inline void* MyCalorimeterHit::operator new(size_t)
{
  return (void *) MyCalorimeterHitAllocator.MallocSingle();
}

inline void MyCalorimeterHit::operator delete(void *aHit)
{
  MyCalorimeterHitAllocator.FreeSingle((MyCalorimeterHit*) aHit);
}

#endif


