//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: GammaRayTelDigi.hh,v 1.1 2001-10-24 13:11:54 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  int operator==(const GammaRayTelDigi&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4int PlaneNumber;    //  (active detector)
  G4int PlaneType;      // (0 or 1 for X or Y plane)
  G4int StripNumber; // strip number
  
public:
  
  inline void SetPlaneNumber(G4int PlaneNum)   {PlaneNumber = PlaneNum;};
  inline void SetPlaneType(G4int PlaneTyp)   {PlaneType = PlaneTyp;};
  inline void SetStripNumber(G4int StripNum)  {StripNumber = StripNum;};
  
  inline G4int GetPlaneNumber() {return PlaneNumber;};
  inline G4int GetPlaneType()   {return PlaneType;};
  inline G4int GetStripNumber() {return StripNumber;};
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<GammaRayTelDigi> GammaRayTelDigitsCollection;

extern G4Allocator<GammaRayTelDigi> GammaRayTelDigiAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* GammaRayTelDigi::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) GammaRayTelDigiAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void GammaRayTelDigi::operator delete(void* aDigi)
{
  GammaRayTelDigiAllocator.FreeSingle((GammaRayTelDigi*) aDigi);
}

#endif









