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

#ifndef G4PiMinusAbsorptionBertini_h
#define G4PiMinusAbsorptionBertini_h 1

// Class Description:
//
// Process for pi- absorption at rest. 
// To be used in your physics list in case you need this physics.

#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4CascadeInterface.hh"
#include "G4HadronicProcessType.hh"

class G4PiMinusAbsorptionBertini : public G4VRestProcess

{ 
private:
    // hide assignment operator as private 
    G4PiMinusAbsorptionBertini& operator=(const G4PiMinusAbsorptionBertini &right);
    G4PiMinusAbsorptionBertini(const G4PiMinusAbsorptionBertini& );
    
public:
    
    G4PiMinusAbsorptionBertini(const G4String& processName ="PiMinusAbsorptionBertini",
                               G4ProcessType   aType = fHadronic);
    
    ~G4PiMinusAbsorptionBertini();
    
    G4bool IsApplicable(const G4ParticleDefinition&);
    
    void PreparePhysicsTable(const G4ParticleDefinition&);
    
    // null physics table
    void BuildPhysicsTable(const G4ParticleDefinition&);
    
    G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 
    
private:
    
    void GenerateSecondaries();
    void PionMinusAbsorption( G4int* );
    
protected:                         
    
    // zero mean lifetime
    G4double GetMeanLifeTime(const G4Track& aTrack,
                             G4ForceCondition* ) 
    {
        G4double result = 0;
        if(aTrack.GetMaterial()->GetNumberOfElements() == 1)
            if(aTrack.GetMaterial()->GetZ()<1.5) result = DBL_MAX;
        return result;
    }
    
    
private:
    // atomic mass of target nucleus
    G4float  currentN;
    
    // charge of target nucleus
    G4float  targetZ;
    
    G4Nucleus targetNucleus;
    G4CascadeInterface* cascade;
    G4ParticleDefinition* pdefPionMinus;
};

#endif

