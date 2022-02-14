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
/// \file InteractionHit.cc
/// \brief Implementation of the CaTS::InteractionHit class

#include <G4VHit.hh>
#include "InteractionHit.hh"
template <class Type>
class G4Allocator;
G4ThreadLocal G4Allocator<InteractionHit>* InteractionHitAllocator = nullptr;

InteractionHit::InteractionHit()
  : G4VHit()
{}

InteractionHit::InteractionHit(G4String pn, G4double p, G4double e, G4double t)
  : G4VHit()
{
  fpname = pn;  // name of secondary particle
  fpmom  = p;   // momentum of secondary particle
  fEkin  = e;   // kinetic energy of secondary particle
  ftheta = t;   // theta of secondary particle
}

InteractionHit::InteractionHit(const InteractionHit& right)
  : G4VHit()
{
  fpname = right.fpname;
  fpmom  = right.fpmom;
  fEkin  = right.fEkin;
  ftheta = right.ftheta;
}

const InteractionHit& InteractionHit::operator=(const InteractionHit& right)
{
  fpname = right.fpname;
  fpmom  = right.fpmom;
  fEkin  = right.fEkin;
  ftheta = right.ftheta;
  return *this;
}

G4int InteractionHit::operator==(const InteractionHit& right) const
{
  return (this == &right) ? 1 : 0;
}
