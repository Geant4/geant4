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
//      ------------ GammaRayTelCalorimeterHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Calorimeter

#ifndef GammaRayTelCalorimeterHit_h
#define GammaRayTelCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelCalorimeterHit : public G4VHit
{
public:
  
  GammaRayTelCalorimeterHit();
  ~GammaRayTelCalorimeterHit();
  GammaRayTelCalorimeterHit(const GammaRayTelCalorimeterHit&);
  const GammaRayTelCalorimeterHit& operator=(const
						GammaRayTelCalorimeterHit&);
  G4bool operator==(const GammaRayTelCalorimeterHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double EdepCAL;  // Energy deposited on the ACD tile
  G4ThreeVector pos; // Position of the hit
  G4int CALBarNumber; // Number of the CAL tile
  G4int CALPlaneNumber;    // Number of the CAL plane
  G4int IsCALPlane;    // Type of the plane (0 X, 1 Y)

public:
  
  inline void AddEnergy(G4double de) {EdepCAL += de;};
  inline void SetCALBarNumber(G4int i) {CALBarNumber = i;};
  inline void SetCALPlaneNumber(G4int i) {CALPlaneNumber = i;};
  inline void SetCALType(G4int i) {IsCALPlane = i;};
  inline void SetPos(G4ThreeVector xyz){ pos = xyz; }
  
  inline G4double GetEdepCAL()     { return EdepCAL; };
  inline G4int    GetCALBarNumber()   { return CALBarNumber; };
  inline G4int    GetCALPlaneNumber()   { return CALPlaneNumber; };
  inline G4int    GetCALType()   {return IsCALPlane;};      
  inline G4ThreeVector GetPos() { return pos; };
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<GammaRayTelCalorimeterHit> GammaRayTelCalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<GammaRayTelCalorimeterHit> *GammaRayTelCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelCalorimeterHit::operator new(size_t)
{
  if (!GammaRayTelCalorimeterHitAllocator)
    GammaRayTelCalorimeterHitAllocator = new G4Allocator<GammaRayTelCalorimeterHit>;
  return (void*) GammaRayTelCalorimeterHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelCalorimeterHit::operator delete(void* aHit)
{
  GammaRayTelCalorimeterHitAllocator->FreeSingle((GammaRayTelCalorimeterHit*) aHit);
}

#endif










