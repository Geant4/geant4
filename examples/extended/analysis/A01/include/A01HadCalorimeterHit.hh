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
// $Id: A01HadCalorimeterHit.hh,v 1.5 2006-06-29 16:31:12 gunter Exp $
// --------------------------------------------------------------
//
#ifndef A01HadCalorimeterHit_h
#define A01HadCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class A01HadCalorimeterHit : public G4VHit
{
  public:

      A01HadCalorimeterHit();
      A01HadCalorimeterHit(G4int iCol,G4int iRow);
      virtual ~A01HadCalorimeterHit();
      A01HadCalorimeterHit(const A01HadCalorimeterHit &right);
      const A01HadCalorimeterHit& operator=(const A01HadCalorimeterHit &right);
      int operator==(const A01HadCalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
      G4int columnID;
      G4int rowID;
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;

  public:
      inline void SetColumnID(G4int z) { columnID = z; }
      inline G4int GetColumnID() const { return columnID; }
      inline void SetRowID(G4int z) { rowID = z; }
      inline G4int GetRowID() const { return rowID; }
      inline void SetEdep(G4double de) { edep = de; }
      inline void AddEdep(G4double de) { edep += de; }
      inline G4double GetEdep() const { return edep; }
      inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
      inline G4ThreeVector GetPos() const { return pos; }
      inline void SetRot(G4RotationMatrix rmat) { rot = rmat; }
      inline G4RotationMatrix GetRot() const { return rot; }

};

typedef G4THitsCollection<A01HadCalorimeterHit> A01HadCalorimeterHitsCollection;

extern G4Allocator<A01HadCalorimeterHit> A01HadCalorimeterHitAllocator;

inline void* A01HadCalorimeterHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)A01HadCalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void A01HadCalorimeterHit::operator delete(void* aHit)
{
  A01HadCalorimeterHitAllocator.FreeSingle((A01HadCalorimeterHit*) aHit);
}

#endif


