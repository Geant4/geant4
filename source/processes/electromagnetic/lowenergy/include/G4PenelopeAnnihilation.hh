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
// Author: Luciano Pandola
//
// History:
// -----------
// 01 Jul 2003   L. Pandola   1st implementation
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics, Positron Annihilation
// Penelope Model
// -------------------------------------------------------------------

#ifndef G4PenelopeAnnihilation_h
#define G4PenelopeAnnihilation_h 1

#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"

class G4PhysicsTableclass;
class G4ParticleDefinition;
class G4VParticleChange;
class G4Track;
class G4Step;
class G4Material;

class G4PenelopeAnnihilation : public G4VRestDiscreteProcess
 {    
 public:  
 
   G4PenelopeAnnihilation(const G4String& processName ="PenAnnih");
 
   ~G4PenelopeAnnihilation();

   G4bool IsApplicable(const G4ParticleDefinition&);
  
   void BuildPhysicsTable(const G4ParticleDefinition&);
  
   void PrintInfoDefinition();

   G4double DumpMeanFreePath(const G4Track& aTrack, 
			     G4double previousStepSize, 
			     G4ForceCondition* condition) 
   { return GetMeanFreePath(aTrack, previousStepSize, condition); }
     
 protected:
   G4double GetMeanFreePath(const G4Track&,G4double,
			    G4ForceCondition*);
     
 public:
 
   G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep); 

   G4double GetMeanLifeTime(const G4Track& aTrack,
			   G4ForceCondition* condition);
              
   G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                                   const G4Step& aStep); 
  
 private:

   G4double calculateCrossSectionPerElectron(G4double energy);

 
 private:
  
   // hide assignment operator as private 
   G4PenelopeAnnihilation& operator=(const G4PenelopeAnnihilation& right);
   G4PenelopeAnnihilation(const G4PenelopeAnnihilation& );
      
 private:
   G4PhysicsTable* meanFreePathTable;

   G4double lowEnergyLimit;      // low  energy limit of the tables
   G4double highEnergyLimit;     // high energy limit of the tables 
   G4int nBins;                  // number of bins in the tables
   G4double cutForLowEnergySecondaryPhotons;
};

#endif
 
