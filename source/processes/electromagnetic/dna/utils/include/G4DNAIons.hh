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
// 31.07.2013: Derived from G4Ions for MT migration (SI)
// Suggested by Makoto et al.

#ifndef G4DNAIons_h
#define G4DNAIons_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                          G4DNAIons                                 ###
// ######################################################################

class G4DNAIons : public G4ParticleDefinition
{
 // Class Description
 //  This is the base class for all Geant4-DNA nuclei, including
 //  charged states 

 protected:
   G4DNAIons(){};


 public:
   G4DNAIons(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable,  G4bool              shortlived,
       const G4String&     subType ="",
       G4int               anti_encoding =0,
       G4double            excitation = 0.0, 
       G4int               isomer = 0
   );

 public:
   virtual    			~G4DNAIons();
   G4DNAIons*    		IonsDefinition();
   G4DNAIons*    		Ions();

 public:
  
  // Get excitation energy of nucleus
  G4double GetExcitationEnergy() const ; 
  
  // Get Isomer level (=0 for ground state)
  G4int GetIsomerLevel() const; 
   
  private:
  G4double theExcitationEnergy; 
  G4int    theIsomerLevel;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
 G4DNAIons* G4DNAIons::IonsDefinition()
{
  return this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
 G4DNAIons* G4DNAIons::Ions() 
{
  return this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
 G4double G4DNAIons::GetExcitationEnergy() const 
{
  return theExcitationEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
 G4int G4DNAIons::GetIsomerLevel() const
{
  return theIsomerLevel;
}
    
#endif








