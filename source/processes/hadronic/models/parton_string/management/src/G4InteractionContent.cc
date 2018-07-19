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
// $Id: G4InteractionContent.cc 100828 2016-11-02 15:25:59Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4InteractionContent----------------
//             by Gunter Folger, June 1998.
//       class for a storing two colliding particles in PartonString Models
// ------------------------------------------------------------

#include "G4InteractionContent.hh"

#include "G4PhysicalConstants.hh" 
#include "G4SystemOfUnits.hh"

G4InteractionContent::G4InteractionContent(G4VSplitableHadron *aPrimaryParticipant)
: theNumberOfHard(0), theNumberOfSoft(0), theNumberOfDiffractive(0),
  theInteractionTime(0.), curStatus(0)
{
  theProjectile=aPrimaryParticipant;
  theTarget=0;
  theProjectileNucleon=0;
  theTargetNucleon=0;
}

G4InteractionContent::~G4InteractionContent()
{}


G4bool G4InteractionContent::operator<(const G4InteractionContent &right) const
{
  return this->GetInteractionTime() < right.GetInteractionTime();
}

void  G4InteractionContent::SetInteractionTime(G4double aValue)
{theInteractionTime = aValue;}

G4double G4InteractionContent::GetInteractionTime() const
{return theInteractionTime;}

void G4InteractionContent::SetStatus(G4int aValue)
{curStatus = aValue;}

G4int G4InteractionContent::GetStatus() const
{return curStatus;}

