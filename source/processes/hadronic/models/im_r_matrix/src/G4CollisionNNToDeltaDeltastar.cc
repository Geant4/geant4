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
// $Id: G4CollisionNNToDeltaDeltastar.cc,v 1.2.2.1 2004/03/24 13:18:33 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToDeltaDeltastar.hh"

#include "G4CollisionNNToDeltaDelta1600.hh"
#include "G4CollisionNNToDeltaDelta1620.hh"
#include "G4CollisionNNToDeltaDelta1700.hh"
#include "G4CollisionNNToDeltaDelta1900.hh"
#include "G4CollisionNNToDeltaDelta1905.hh"
#include "G4CollisionNNToDeltaDelta1910.hh"
#include "G4CollisionNNToDeltaDelta1920.hh"
#include "G4CollisionNNToDeltaDelta1930.hh"
#include "G4CollisionNNToDeltaDelta1950.hh"

typedef
GROUP9(G4CollisionNNToDeltaDelta1600, 
      G4CollisionNNToDeltaDelta1620,
      G4CollisionNNToDeltaDelta1700,
      G4CollisionNNToDeltaDelta1900,
      G4CollisionNNToDeltaDelta1905,
      G4CollisionNNToDeltaDelta1910,
      G4CollisionNNToDeltaDelta1920,
      G4CollisionNNToDeltaDelta1930,
      G4CollisionNNToDeltaDelta1950) theChannels;
      
G4CollisionNNToDeltaDeltastar::G4CollisionNNToDeltaDeltastar()
{ 
  Register aR;
  G4ForEach<theChannels>::Apply(&aR, this); 
}

