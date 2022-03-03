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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file DRCalorimeterHit.cc
/// \brief Implementation of the CaTS::DRCalorimeterHit class

// Geant4 headers
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
// project headers
#include "DRCalorimeterHit.hh"
G4ThreadLocal G4Allocator<DRCalorimeterHit>* DRCalorimeterHitAllocator =
  nullptr;

DRCalorimeterHit::DRCalorimeterHit()
  : G4VHit()
{}

DRCalorimeterHit::DRCalorimeterHit(unsigned int i, G4double e, G4double em,
                                   unsigned int nc, G4double t, G4ThreeVector p)
  : G4VHit()
{
  fid       = i;
  fEdep     = e;
  fem_Edep  = em;
  fNceren   = nc;
  ftime     = t;
  fposition = p;
}

DRCalorimeterHit::DRCalorimeterHit(const DRCalorimeterHit& right)
  : G4VHit()
{
  this->fid       = right.fid;
  this->fEdep     = right.fEdep;
  this->fem_Edep  = right.fem_Edep;
  this->fNceren   = right.fNceren;
  this->ftime     = right.ftime;
  this->fposition = right.fposition;
}

const DRCalorimeterHit& DRCalorimeterHit::operator=(
  const DRCalorimeterHit& right)
{
  this->fid       = right.fid;
  this->fEdep     = right.fEdep;
  this->fem_Edep  = right.fem_Edep;
  this->fNceren   = right.fNceren;
  this->ftime     = right.ftime;
  this->fposition = right.fposition;
  return *this;
}

G4bool DRCalorimeterHit::operator==(const DRCalorimeterHit& right) const
{
  return (this == &right) ? true : false;
}

void DRCalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fposition);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1., 0., 0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}
