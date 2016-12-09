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
// $Id: G4LPhysicsFreeVector.cc 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
// --------------------------------------------------------------------
// Class G4LPhysicsFreeVector
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data
//
// F.W. Jones, TRIUMF, 04-JUN-96
// 
// Modified:
//    19 Jun. 2009, V.Ivanchenko : removed hidden bin 
//    19 Jun. 2009, V.Ivanchenko : removed FindBinLocation 
//
// --------------------------------------------------------------------

#include "G4LPhysicsFreeVector.hh"

// --------------------------------------------------------------------

G4LPhysicsFreeVector::G4LPhysicsFreeVector()
   : G4PhysicsFreeVector()
{}

// --------------------------------------------------------------------

G4LPhysicsFreeVector::G4LPhysicsFreeVector(size_t length, G4double, G4double)
  : G4PhysicsFreeVector(length)
{}  

// --------------------------------------------------------------------

G4LPhysicsFreeVector::~G4LPhysicsFreeVector()
{}

// --------------------------------------------------------------------
