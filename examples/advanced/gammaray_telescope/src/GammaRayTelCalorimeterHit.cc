// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelCalorimeterHit.cc,v 1.1 2001-03-05 13:58:22 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelCalorimeterHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelCalorimeterHit.hh"

G4Allocator<GammaRayTelCalorimeterHit> GammaRayTelCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterHit::GammaRayTelCalorimeterHit()
{
  EdepCAL = 0.; 
  CALBarNumber = 0;
  CALPlaneNumber = 0;
  IsCALPlane = 0;
  pos = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterHit::~GammaRayTelCalorimeterHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterHit::GammaRayTelCalorimeterHit(const GammaRayTelCalorimeterHit& right)
{
  EdepCAL = right.EdepCAL; 
  CALBarNumber = right.CALBarNumber;
  CALPlaneNumber = right.CALPlaneNumber;
  IsCALPlane = right.IsCALPlane;
  pos = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const GammaRayTelCalorimeterHit& GammaRayTelCalorimeterHit::operator=(const GammaRayTelCalorimeterHit& right)
{
  EdepCAL = right.EdepCAL; 
  CALBarNumber = right.CALBarNumber;
  CALPlaneNumber = right.CALPlaneNumber;
  IsCALPlane = right.IsCALPlane;
  pos = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int GammaRayTelCalorimeterHit::operator==(const GammaRayTelCalorimeterHit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











