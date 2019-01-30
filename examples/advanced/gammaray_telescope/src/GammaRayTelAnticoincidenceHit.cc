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
//      ------------ GammaRayTelAnticoincidenceHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelAnticoincidenceHit.hh"

G4ThreadLocal G4Allocator<GammaRayTelAnticoincidenceHit> *GammaRayTelAnticoincidenceHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceHit::GammaRayTelAnticoincidenceHit()
{
  EdepACD = 0.; 
  ACDTileNumber = 0; 
  IsACDPlane = 0;
  pos = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceHit::~GammaRayTelAnticoincidenceHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceHit::GammaRayTelAnticoincidenceHit(const GammaRayTelAnticoincidenceHit& right)
  : G4VHit()
{
  EdepACD = right.EdepACD; 
  ACDTileNumber = right.ACDTileNumber;
  IsACDPlane = right.IsACDPlane;
  pos = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const GammaRayTelAnticoincidenceHit& GammaRayTelAnticoincidenceHit::operator=(const GammaRayTelAnticoincidenceHit& right)
{
  EdepACD = right.EdepACD; 
  ACDTileNumber = right.ACDTileNumber;
  IsACDPlane = right.IsACDPlane;
  pos = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelAnticoincidenceHit::operator==(const GammaRayTelAnticoincidenceHit& right) const
{
   return((EdepACD==right.EdepACD)&&(ACDTileNumber==right.ACDTileNumber)&&(IsACDPlane==right.IsACDPlane)&& (pos==right.pos));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceHit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceHit::Print()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











