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

#include "globals.hh"
#include "G4CollisionNNToDeltaDelta1950.hh"
#include "G4ConcreteNNToDeltaDeltastar.hh"

G4CollisionNNToDeltaDelta1950::G4CollisionNNToDeltaDelta1950()
{ 
	components=0;
	MakeNNToDeltaDelta<Dm_1950PC, D0_1950PC, Dp_1950PC, Dpp_1950PC,
	                   G4ConcreteNNToDeltaDeltastar>::Make(this);
}

G4CollisionNNToDeltaDelta1950::~G4CollisionNNToDeltaDelta1950()
{
	if (components) delete components;
	components=0;
}



