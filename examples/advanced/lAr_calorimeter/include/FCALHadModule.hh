//----------------------------------------------------------------------------
//
//   Name of file:      FCALHadModule.hh
//   Author:            Mathieu Fontaine           Rachid Mazini
//                      fontaine@lps.umontreal.ca  Rachid.Mazini@cern.ch
//   Language:          C++
//   Tested on:         g++
//   Prerequisites:     None
//   Purpose:           Header file for FCALHadModule.cc, which defines
//                      the  geometry of the FCAL HadModule 0.
//   Developped:        10-March-2000   M.F.
//
//----------------------------------------------------------------------------

#ifndef FCALHadModule_h
#define FCALHadModule_h 1

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "FCALHadModuleSD.hh"

class FCALHadModule
{
public:

  FCALHadModule();
  ~FCALHadModule();

public:

  G4LogicalVolume * Construct();

  void InitializeGeometry(); 
  G4int GetF2TileID(G4int); 

private:
  G4int NF2LarGap;
  G4int* F2LArGapID;
  G4int* F2LArIX;
  G4int* F2LArJY;
  G4int* F2LArITile;
  G4double* F2LArGapPosX;
  G4double* F2LArGapPosY;

private:

  G4double HadModuleRMin, HadModuleRMax, HadModuleLenght;
  G4double HadModuleStartPhi, HadModuleDPhi;
  G4double FCAL2HadSmart;

  G4double WAbsorberRMin, WAbsorberRMax, WAbsorberLenght;
  G4double WAbsorberStartPhi, WAbsorberDPhi;

  G4double CuPlateLenght, CuPlateAPosZ, CuPlateBPosZ;

  G4int NCableTroff;
  G4double F2TroffRmin, F2TroffRmax, F2TroffMainLenght;
  G4double F2TroffABLenght, F2TroffStartPhi,  F2TroffDphi; 
  G4double F2TroffRotZ, F2TroffABPosZ; 

  G4double F2LArGapRmin, F2LArGapRmax, F2LArGapLenght;
  G4double F2LArGapStartPhi, F2LArGapDphi;
  G4double F2RodRmin, F2RodRmax, F2RodLenght, F2RodStartPhi, F2RodDphi;

  FCALHadModuleSD* FcalHadModuleSD;

};

#endif  /* FCALHadModule.hh */
