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
#ifndef G4FieldPropagation_h
#define G4FieldPropagation_h 1

#include "G4KineticTrackVector.hh"

class G4FieldPropagation 
{
public:
   G4FieldPropagation() {}
   G4FieldPropagation(const G4FieldPropagation &) {}

   virtual ~G4FieldPropagation() {}

   // Operators
   const G4FieldPropagation & operator=(const G4FieldPropagation &right);

   int operator==(const G4FieldPropagation &right) const;
   int operator!=(const G4FieldPropagation &right) const;

    // Methods

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


