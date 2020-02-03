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
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelAnticoincidenceHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Anticoincidence

#ifndef GammaRayTelAnticoincidenceHit_h
#define GammaRayTelAnticoincidenceHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelAnticoincidenceHit : public G4VHit
{
public:
  
  GammaRayTelAnticoincidenceHit();
  ~GammaRayTelAnticoincidenceHit();
  GammaRayTelAnticoincidenceHit(const GammaRayTelAnticoincidenceHit&);
  const GammaRayTelAnticoincidenceHit& operator=(const
						GammaRayTelAnticoincidenceHit&);
  G4bool operator==(const GammaRayTelAnticoincidenceHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double EdepACD;  // Energy deposited on the ACD tile
  G4ThreeVector pos; // Position of the hit
  G4int ACDTileNumber; // Number of the ACD tile
  G4int IsACDPlane;    // Type of the plane (0 top, 1 L-R, 2 F-R)

public:
  
  inline void AddEnergy(G4double de) {EdepACD += de;};
  inline void SetACDTileNumber(G4int i) {ACDTileNumber = i;};
  inline void SetACDType(G4int i) {IsACDPlane = i;};
  inline void SetPos(G4ThreeVector xyz){ pos = xyz; }
  
  inline G4double GetEdepACD()     { return EdepACD; };
  inline G4int    GetACDTileNumber()   { return ACDTileNumber; };
  inline G4int    GetACDType()   {return IsACDPlane;};      
  inline G4ThreeVector GetPos() { return pos; };
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<GammaRayTelAnticoincidenceHit> GammaRayTelAnticoincidenceHitsCollection;

extern G4ThreadLocal G4Allocator<GammaRayTelAnticoincidenceHit> *GammaRayTelAnticoincidenceHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelAnticoincidenceHit::operator new(size_t)
{
  if (!GammaRayTelAnticoincidenceHitAllocator)
    GammaRayTelAnticoincidenceHitAllocator = new G4Allocator<GammaRayTelAnticoincidenceHit> ;
  return (void*) GammaRayTelAnticoincidenceHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelAnticoincidenceHit::operator delete(void* aHit)
{
  GammaRayTelAnticoincidenceHitAllocator->FreeSingle((GammaRayTelAnticoincidenceHit*) aHit);
}

#endif










