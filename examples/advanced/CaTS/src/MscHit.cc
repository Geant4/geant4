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
/// \file MscHit.cc
/// \brief Implementation of the CaTS::MscHit class

// Geant4 headers
#include <G4ThreeVector.hh>
#include <G4VHit.hh>
// project headers
#include "MscHit.hh"

template <class Type>
class G4Allocator;
G4ThreadLocal G4Allocator<MscHit>* MscHitAllocator = nullptr;

MscHit::MscHit()
  : G4VHit()
{}
MscHit::MscHit(G4double E, G4ThreeVector mom)
  : G4VHit()
{
  this->fkinE     = E;
  this->fmomentum = mom;
}

MscHit::MscHit(const MscHit& right)
  : G4VHit()
{
  fkinE     = right.fkinE;
  fmomentum = right.fmomentum;
}

const MscHit& MscHit::operator=(const MscHit& right)
{
  fkinE     = right.fkinE;
  fmomentum = right.fmomentum;
  return *this;
}

G4bool MscHit::operator==(const MscHit& right) const
{
  return (this == &right) ? true : false;
}

void MscHit::Draw() {}
