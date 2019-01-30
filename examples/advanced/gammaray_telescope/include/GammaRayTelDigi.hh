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
//      ------------ GammaRayTelDigi  ------
//           by F.Longo, R.Giannitrapani & G.Santin (24 oct 2001)
//
// ************************************************************
// This Class describe the digits 

#ifndef GammaRayTelDigi_h
#define GammaRayTelDigi_h 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelDigi : public G4VDigi
{

public:
  
  GammaRayTelDigi();
  ~GammaRayTelDigi();
  GammaRayTelDigi(const GammaRayTelDigi&);
  const GammaRayTelDigi& operator=(const GammaRayTelDigi&);
  G4bool operator==(const GammaRayTelDigi&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4int PlaneNumber;    //  (active detector)
  G4int PlaneType;      // (0 or 1 for X or Y plane)
  G4int StripNumber; // strip number
  G4int DigiType;        // (0 == TKR, 1 == CAL, 2 == ACD)
  G4double Energy; // only for CAL 
  
public:
  
  inline void SetPlaneNumber(G4int PlaneNum)   {PlaneNumber = PlaneNum;};
  inline void SetPlaneType(G4int PlaneTyp)   {PlaneType = PlaneTyp;};
  inline void SetStripNumber(G4int StripNum)  {StripNumber = StripNum;};
  inline void SetDigiType(G4int DigiID)  {DigiType = DigiID;};
  inline void SetEnergy(G4double Ene)  {Energy = Ene;};

  inline G4int GetPlaneNumber() {return PlaneNumber;};
  inline G4int GetPlaneType()   {return PlaneType;};
  inline G4int GetStripNumber() {return StripNumber;};
  inline G4int GetDigiType() {return DigiType;};
  inline G4double GetEnergy()  {return Energy;};
  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<GammaRayTelDigi> GammaRayTelDigitsCollection;

extern G4ThreadLocal G4Allocator<GammaRayTelDigi> *GammaRayTelDigiAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelDigi::operator new(size_t)
{
  if (!GammaRayTelDigiAllocator)
    GammaRayTelDigiAllocator = new G4Allocator<GammaRayTelDigi>;
  return (void*) GammaRayTelDigiAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelDigi::operator delete(void* aDigi)
{
  GammaRayTelDigiAllocator->FreeSingle((GammaRayTelDigi*) aDigi);
}

#endif









