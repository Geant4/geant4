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
/// \file ExN04CalorimeterHit.hh
/// \brief Definition of the ExN04CalorimeterHit class
//

#ifndef ExN04CalorimeterHit_h
#define ExN04CalorimeterHit_h 1


//MSH_include_begin
#include "MarshaledG4String.h"
//MSH_include_end


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"


//MSH_BEGIN
class ExN04CalorimeterHit : public G4VHit
{
  public:

      ExN04CalorimeterHit();
      ExN04CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi);
      ~ExN04CalorimeterHit();
      ExN04CalorimeterHit(const ExN04CalorimeterHit &right);
      const ExN04CalorimeterHit& operator=(const ExN04CalorimeterHit &right);
      G4bool operator==(const ExN04CalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:

      G4int ZCellID; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetZ(); }]
  [elementSet: { $THIS->SetCellID($ELEMENT, $THIS->GetPhi()); }] */

      G4int PhiCellID; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetPhi(); }]
  [elementSet: { $THIS->SetCellID($THIS->GetZ(), $ELEMENT); }] */ 


    G4double edep; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetEdep(); }]
  [elementSet: { $THIS->SetEdep($ELEMENT); }] */

      G4ThreeVector pos; /*MSH: primitive
  [elementGet: { $ELEMENT = $THIS->GetPos(); }]
  [elementSet: { $THIS->SetPos($ELEMENT); }] */

  G4RotationMatrix rot; 
  const G4LogicalVolume* pLogV;



  public:
      inline void SetCellID(G4int z,G4int phi)
      {
        ZCellID = z;
        PhiCellID = phi;
      }
      inline G4int GetZ() { return ZCellID; }
      inline G4int GetPhi() { return PhiCellID; }
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de)
      { edep += de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetRot(G4RotationMatrix rmat)
      { rot = rmat; }
      inline G4RotationMatrix GetRot()
      { return rot; }
      inline const G4LogicalVolume * GetLogV()
      { return pLogV; }

};
//MSH_END

typedef G4THitsCollection<ExN04CalorimeterHit> ExN04CalorimeterHitsCollection;

extern G4Allocator<ExN04CalorimeterHit> ExN04CalorimeterHitAllocator;

inline void* ExN04CalorimeterHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExN04CalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04CalorimeterHit::operator delete(void *aHit)
{
  ExN04CalorimeterHitAllocator.FreeSingle((ExN04CalorimeterHit*) aHit);
}

#endif


