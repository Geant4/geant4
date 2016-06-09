//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

