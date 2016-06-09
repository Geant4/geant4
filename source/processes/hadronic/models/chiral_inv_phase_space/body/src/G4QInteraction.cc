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
// $Id: G4QInteraction.cc,v 1.2 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// ------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QInteraction----------------
//            Created by Mikhail Kossov Oct, 2006
//   class for a storing colliding particles in PartonString Models
//   For comparison mirror member functions are taken from G4 class:
//   G4InteractionContent
// -------------------------------------------------------------------

#include "G4QInteraction.hh"

G4QInteraction::G4QInteraction(G4QHadron* aProjectile) : theProjectile(aProjectile),
  theTarget(0),theNumberOfHard(0),theNumberOfSoft(0),theNumberOfDiffractive(0)
{}

G4QInteraction::G4QInteraction(const G4QInteraction &right) :
  theProjectile(right.GetProjectile()), theTarget(right.GetTarget()),
  theNumberOfHard(0), theNumberOfSoft(0), theNumberOfDiffractive(0)
{}

G4QInteraction::~G4QInteraction()
{}


