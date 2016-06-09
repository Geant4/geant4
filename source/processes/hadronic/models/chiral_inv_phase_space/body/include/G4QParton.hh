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

#ifndef G4QParton_h
#define G4QParton_h 1

// $Id$
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4QParton ----------------
//      by Mikhail Kosov, October 2006
//      class for Quark-Parton for quark-level models
//   For comparison mirror member functions are taken from G4 class:
//   G4Parton
// -----------------------------------------------------------------
// Short description: The Quark-Gluon String consists of the partons, which
// are quarks and some times gluons.
// ------------------------------------------------------------------------

#include "globals.hh"
#include "G4QContent.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include <iostream>
#include "Randomize.hh"

class G4QParton
{
 public:
  // Constructors
  G4QParton();                                            // Default fullRandom constructor
  G4QParton(G4int aPGG);                                  // Collor/Spin are still random
  G4QParton(const G4QParton &right);
  G4QParton(const G4QParton* right);

  // Destructor
  ~G4QParton();

  // Operators
  const G4QParton& operator=(const G4QParton &right);
  G4int operator==(const G4QParton &right) const      {return this==&right;} 
  G4int operator!=(const G4QParton &right) const      {return this!=&right;} 

  // Modifiers
  void DefineEPz(G4LorentzVector hadr4Mom){theMomentum+=hadr4Mom*theX;} // CHIPS solution
  void DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection);
  void SetPDGCode(G4int aPDG);
  void SetColour(G4int aColour)                       {theColour = aColour;}
  void SetX(G4double anX)                             {theX = anX;}
  void Set4Momentum(const G4LorentzVector& aMomentum) {theMomentum=aMomentum;}
  void SetPosition(const G4ThreeVector& aPosition)    {thePosition=aPosition;}
  void SetSpinZ(G4double aSpinZ)                      {theSpinZ = aSpinZ;}
  G4bool ReduceDiQADiQ(G4QParton* d1, G4QParton* d2); // Reduce DiQ-aDiQ to Q-aQ (general)

  // Selectors
  G4int GetPDGCode() const                            {return PGGCode;}
  G4QContent GetQC() const                            {return QCont;}
  const G4ThreeVector& GetPosition() const            {return thePosition;}
  const G4LorentzVector& Get4Momentum() const         {return theMomentum;} 
  G4double GetX()                                     {return theX;}    
  G4int GetColour()                                   {return theColour;}    
  G4double GetSpinZ()                                 {return theSpinZ;}
  const G4int& GetType() const                        {return theType;}
 private: 
  // Body 
  G4int                 PGGCode;
  G4QContent            QCont;   // Quark Content of the parton
  G4int                 theType; // 0 = gluon, 1 = quark-antiquark, 2 = diquark/antidiquark
  G4int                 theColour;
  G4double              theSpinZ;
  G4double              theX;
  G4ThreeVector         thePosition;
  G4LorentzVector       theMomentum;
};

#endif
