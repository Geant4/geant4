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
#ifndef G4FieldPropagation_h
#define G4FieldPropagation_h 1

#include "G4KineticTrackVector.hh"

class G4FieldPropagation 
{
public:
   G4FieldPropagation() {}
   G4FieldPropagation(const G4FieldPropagation &) {}

   virtual ~G4FieldPropagation();

private:   // Operators
   const G4FieldPropagation & operator=(const G4FieldPropagation &right);

   int operator==(const G4FieldPropagation &right) const;
   int operator!=(const G4FieldPropagation &right) const;

public:    // Methods

   // only theActive are propagated, nothing else
   // only theSpectators define the field, nothing else
   virtual void Transport(G4KineticTrackVector &theActive, const G4KineticTrackVector &theSpectators, G4double theTimeStep) = 0;

   virtual G4double GetExcitationEnergy(G4int nHit, const G4KineticTrackVector &theParticles) = 0;
   
   // methods for calculating potentials for different types of particles
   virtual void Init(G4int z, G4int a) = 0; // prepare potentials' functions
   
   // aPosition is relative to the nucleus center
   virtual G4double GetNeutronPotential(G4double radius) = 0;
   virtual G4double GetNeutronPotential(G4ThreeVector &aPosition) = 0;
   
   virtual G4double GetProtonPotential(G4double radius) = 0;
   virtual G4double GetProtonPotential(G4ThreeVector &aPosition) = 0;
    
   virtual G4double GetAntiprotonPotential(G4double radius) = 0;
   virtual G4double GetAntiprotonPotential(G4ThreeVector &aPosition) = 0;

   virtual G4double GetKaonPotential(G4double radius) = 0;
   virtual G4double GetKaonPotential(G4ThreeVector &aPosition) = 0;
   
   virtual G4double GetPionPotential(G4double radius) = 0;
   virtual G4double GetPionPotential(G4ThreeVector &aPosition) = 0;
};

#endif // G4FieldPropagation_h


