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
// $Id: G4AlphaCoulombBarrier.cc 97097 2016-05-25 07:50:31Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4AlphaCoulombBarrier.hh"

G4AlphaCoulombBarrier::G4AlphaCoulombBarrier() : G4CoulombBarrier(4,2) {}
G4AlphaCoulombBarrier::~G4AlphaCoulombBarrier() {}

G4double G4AlphaCoulombBarrier::BarrierPenetrationFactor(G4int aZ) const
{
  // Data comes from 
  // Dostrovsky, Fraenkel and Friedlander
  // Physical Review, vol 116, num. 3 1959
  // 
  // const G4int size = 5;
  // const G4double Zlist[size]  = {10.0, 20.0, 30.0, 50.0, 70.0};
  // const G4double Kalpha[size] = {0.68, 0.82, 0.91, 0.97, 0.98};
  return  (aZ >= 70) ? 0.98 :
    (((0.23684e-5*aZ) - 0.42143e-3)*aZ + 0.25222e-1)*aZ + 0.46699;
}
