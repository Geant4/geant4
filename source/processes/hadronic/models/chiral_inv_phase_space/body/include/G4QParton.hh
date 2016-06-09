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

// $Id: G4QParton.hh,v 1.2 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-01 $
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

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include <iostream>
#include "G4ParticleTable.hh"
#include "Randomize.hh"

class G4QParton
{
 public:
  // Constructors
  G4QParton()
  {
    // CHIPS is working only with u, d, and s quarks (SU(3)xSU(3)) (no gluons! M.K.)
    // Random Flavor/Colour/Spin definition for default constructor (with .3 s-suppresion)
    PDGencoding=(G4int)(2.3*G4UniformRand())+1; //@@ What about antiquarks? (M.K.)
   	theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);
   	// random colour (1,2,3)=(R,G,B) for quarks and (-1,-2,-3)=(aR,aG,aB) for anti-quarks
    theColour = (G4int)(3*G4UniformRand())+1;
    if(theColour>3) theColour = 3;
    theIsoSpinZ = theDefinition->GetPDGIsospin3();
    theSpinZ = (G4int)(2*G4UniformRand()) - 0.5;
  }
		G4QParton(G4int PDGencoding);                            // Collor/Spin are still random
  G4QParton(const G4QParton &right);
  G4QParton(const G4QParton* right);

  // Destructor
  ~G4QParton();

  // Operators
  const G4QParton& operator=(const G4QParton &right);
  G4int operator==(const G4QParton &right) const      {return this==&right;}	
  G4int operator!=(const G4QParton &right) const      {return this!=&right;}	

  // Modifiers
  void DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection);
  void SetColour(G4int aColour)                       {theColour = aColour;}
  void SetX(G4double anX)                             {theX = anX;}
  void Set4Momentum(const G4LorentzVector& aMomentum) {theMomentum=aMomentum;}
  void SetPosition(const G4ThreeVector& aPosition)    {thePosition=aPosition;}
  void SetIsoSpinZ(G4double anIsoSpinZ)               {theIsoSpinZ = anIsoSpinZ;}
  void SetSpinZ(G4double aSpinZ)                      {theSpinZ = aSpinZ;}

  // Selectors
  G4int GetPDGCode() const                            {return PDGencoding;}
  G4ParticleDefinition* GetDefinition()               {return theDefinition;}    
  const G4ThreeVector& GetPosition()const             {return thePosition;}
  const G4LorentzVector& Get4Momentum() const         {return theMomentum;} 
  G4double GetX()                                     {return theX;}    
  G4int GetColour()                                   {return theColour;}    
  G4double GetSpinZ()                                 {return theSpinZ;}
  G4double GetIsoSpinZ()                              {return theIsoSpinZ;}
  G4double GetMass()                                  {return theDefinition->GetPDGMass();}
  G4String GetParticleSubType()               {return theDefinition->GetParticleSubType();}
 private: 
		// Body 
  G4int                 PDGencoding;
  G4ParticleDefinition* theDefinition;
  G4LorentzVector       theMomentum;
  G4ThreeVector         thePosition;
  G4int                 theColour;
  G4double              theIsoSpinZ;
  G4double              theSpinZ;
  G4double              theX;   
};

#endif
