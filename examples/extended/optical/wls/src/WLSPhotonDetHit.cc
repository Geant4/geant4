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
//

#include "WLSPhotonDetHit.hh"

G4Allocator<WLSPhotonDetHit> WLSPhotonDetHitAllocator;

WLSPhotonDetHit::WLSPhotonDetHit()
{
  posExit     = G4ThreeVector(0., 0., 0.);
  posArrive   = G4ThreeVector(0., 0., 0.);
  arrivalTime = 0.;
}

WLSPhotonDetHit::WLSPhotonDetHit(G4ThreeVector pExit,
                                 G4ThreeVector pArrive,
                                 G4double pTime)
{
  posExit     = pExit;
  posArrive   = pArrive;
  arrivalTime = pTime;
}

WLSPhotonDetHit::~WLSPhotonDetHit() { }

WLSPhotonDetHit::WLSPhotonDetHit(const WLSPhotonDetHit &right)
  : G4VHit()
{
  *this = right;
}

const WLSPhotonDetHit& WLSPhotonDetHit::operator=(const WLSPhotonDetHit &right)
{
  posExit     = right.posExit;
  posArrive   = right.posArrive;
  arrivalTime = right.arrivalTime;

  return *this;
}

G4int WLSPhotonDetHit::operator==(const WLSPhotonDetHit& right) const
{
  return posExit     == right.posExit    &&
         posArrive   == right.posArrive  &&
	 arrivalTime == right.arrivalTime;  
}

void WLSPhotonDetHit::Draw(){ }

void WLSPhotonDetHit::Print(){ }
