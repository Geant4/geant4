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

#include "G4CollisionNNToNDeltastar.hh"

#include "globals.hh"
#include "G4ConcreteNNToNDeltaStar.hh"
#include "G4CollisionNNToNDelta1600.hh"
#include "G4CollisionNNToNDelta1620.hh"
#include "G4CollisionNNToNDelta1700.hh"
#include "G4CollisionNNToNDelta1900.hh"
#include "G4CollisionNNToNDelta1905.hh"
#include "G4CollisionNNToNDelta1910.hh"
#include "G4CollisionNNToNDelta1920.hh"
#include "G4CollisionNNToNDelta1930.hh"
#include "G4CollisionNNToNDelta1950.hh"

typedef
GROUP9(G4CollisionNNToNDelta1600, 
      G4CollisionNNToNDelta1620,
      G4CollisionNNToNDelta1700,
      G4CollisionNNToNDelta1900,
      G4CollisionNNToNDelta1905,
      G4CollisionNNToNDelta1910,
      G4CollisionNNToNDelta1920,
      G4CollisionNNToNDelta1930,
      G4CollisionNNToNDelta1950) theChannels;
      
G4CollisionNNToNDeltastar::G4CollisionNNToNDeltastar()
{ 
  Register aR;
  G4ForEach<theChannels>::Apply(&aR, this); 
}

