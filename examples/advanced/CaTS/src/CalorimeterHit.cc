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
/// \file CalorimeterHit.cc
/// \brief Implementation of the CaTS::CalorimeterHit class

// Geant4 headers
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
// project headers
#include "CalorimeterHit.hh"
#include "G4VVisManager.hh"
G4ThreadLocal G4Allocator<CalorimeterHit>* CalorimeterHitAllocator = nullptr;

CalorimeterHit::CalorimeterHit()
  : G4VHit()
{}

CalorimeterHit::CalorimeterHit(unsigned i, G4double e, G4double em, G4double t,
                               G4ThreeVector p)
  : G4VHit()
{
  fid       = i;
  fEdep     = e;
  fem_Edep  = em;
  ftime     = t;
  fposition = p;
}

CalorimeterHit::CalorimeterHit(const CalorimeterHit& right)
  : G4VHit()
{
  // No need to call class member functions via this->
  fid       = right.fid;
  fEdep     = right.fEdep;
  fem_Edep  = right.fem_Edep;
  ftime     = right.ftime;
  fposition = right.fposition;
}

const CalorimeterHit& CalorimeterHit::operator=(const CalorimeterHit& right)
{
  // No need to call class member functions via this->
  fid       = right.fid;
  fEdep     = right.fEdep;
  fem_Edep  = right.fem_Edep;
  ftime     = right.ftime;
  fposition = right.fposition;
  return *this;
}

G4bool CalorimeterHit::operator==(const CalorimeterHit& right) const
{
  return (this == &right) ? true : false;
}

void CalorimeterHit::Draw()
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
