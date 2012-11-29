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

#ifndef G4QInteraction_h
#define G4QInteraction_h 1

// -------------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4QInteraction----------------
//                 Created by Mikhail Kossov Oct,2006
//   class for a storing colliding particles in PartonString Models
//   For comparison mirror member functions are taken from G4 class:
//   G4InteractionContent
// -------------------------------------------------------------------
//  Short description: Classify the interaction in soft/hard/diffractive
//  parts for firther treatment by the QGS algorithm. Pure data class...
//  Except for the SplitHadrons, which can be done elsewhere (@@ M.K.)
// ---------------------------------------------------------------------

#include "globals.hh"
#include "G4QHadron.hh"

class G4QInteraction 
{
public:
  G4QInteraction(G4QHadron* aPrimaryParticipant);
  ~G4QInteraction();
  int operator==(const G4QInteraction &right) const     {return this==&right;}
  int operator!=(const G4QInteraction &right) const     {return this!=&right;}
  G4QHadron* GetProjectile() const                      {return theProjectile;}
  G4QHadron* GetTarget() const                          {return theTarget;}
  G4int GetNumberOfSoftCollisions()                     {return theNumberOfSoft;}
  G4int GetNumberOfDINRCollisions()                     {return theNumberOfDINR;}
  G4int GetNumberOfHardCollisions()                     {return theNumberOfHard;}
  G4int GetNumberOfDiffractiveCollisions()              {return theNumberOfDiffractive;}
  void  SetTarget(G4QHadron* aTarget)                   {theTarget = aTarget;}
  void  SetProjectile(G4QHadron* aProjectile)           {theProjectile = aProjectile;}
  void  SetNumberOfSoftCollisions(G4int nofSoft)        {theNumberOfSoft = nofSoft;}
  void  SetNumberOfDINRCollisions(G4int nofDINR)        {theNumberOfDINR = nofDINR;}
  void  SetNumberOfHardCollisions(G4int nofHard)        {theNumberOfHard = nofHard;}
  void  SetNumberOfDiffractiveCollisions(G4int nofDiff) {theNumberOfDiffractive = nofDiff;}
  void  SplitHadrons()
            {if(theProjectile)theProjectile->SplitUp(); if(theTarget)theTarget->SplitUp();}
public:
  G4QInteraction(){}
  G4QInteraction(const G4QInteraction &right);
  const G4QInteraction & operator=(const G4QInteraction &right);
private:
  // Body
  G4QHadron* theProjectile;
  G4QHadron* theTarget;
  G4int theNumberOfDINR;
  G4int theNumberOfHard;
  G4int theNumberOfSoft;
  G4int theNumberOfDiffractive;
};

#endif
