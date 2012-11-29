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
// $Id$
//
// ------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QInteraction----------------
//            Created by Mikhail Kossov Oct, 2006
//   class for colliding particles (hadrons) in Parton String Models
//   For comparison mirror member function is G4InteractionContent
// ---------------------------------------------------------------------
//  Short description: Classify the interaction in soft/hard/diffractive
//  parts for firther treatment by the QGS algorithm.
// ---------------------------------------------------------------------

#include "G4QInteraction.hh"

G4QInteraction::G4QInteraction(G4QHadron* aProjectile) :
  theProjectile(aProjectile), theTarget(0), theNumberOfDINR(0),
  theNumberOfHard(0),theNumberOfSoft(0),theNumberOfDiffractive(0)
{}

G4QInteraction::G4QInteraction(const G4QInteraction &right) :
  theProjectile(right.GetProjectile()), theTarget(right.GetTarget()),
  theNumberOfDINR(0), theNumberOfHard(0), theNumberOfSoft(0), theNumberOfDiffractive(0)
{}

G4QInteraction::~G4QInteraction()
{
  //delete theProjectile;
  //if(theTarget) delete theTarget;
}
