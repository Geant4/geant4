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
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelCalorimeterHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelCalorimeterHit.hh"

G4ThreadLocal G4Allocator<GammaRayTelCalorimeterHit> *GammaRayTelCalorimeterHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterHit::GammaRayTelCalorimeterHit()
{
  EdepCAL = 0.; 
  CALBarNumber = 0;
  CALPlaneNumber = 0;
  IsCALPlane = 0;
  pos = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterHit::~GammaRayTelCalorimeterHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelCalorimeterHit::GammaRayTelCalorimeterHit(const GammaRayTelCalorimeterHit& right)
  :G4VHit()
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

G4bool GammaRayTelCalorimeterHit::operator==(const GammaRayTelCalorimeterHit& right) const
{
  return((EdepCAL==right.EdepCAL)&&(CALBarNumber==right.CALBarNumber)&&(CALPlaneNumber==right.CALPlaneNumber)&&(IsCALPlane==right.IsCALPlane)&& (pos==right.pos));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelCalorimeterHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











