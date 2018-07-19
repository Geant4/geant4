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
// $Id: G4EnergyRangeManager.hh 98067 2016-07-01 16:33:54Z gcosmo $
//
 // Hadronic Process: Energy Range Manager
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Last modified: 24-Mar-1997
 
#ifndef G4EnergyRangeManager_h
#define G4EnergyRangeManager_h 1
 
#include "G4HadronicInteraction.hh"
#include <vector>
 
class G4EnergyRangeManager 
{     
public:
    
  explicit G4EnergyRangeManager();
 
  ~G4EnergyRangeManager();
    
  G4EnergyRangeManager(const G4EnergyRangeManager& right);
    
  G4EnergyRangeManager& operator=( const G4EnergyRangeManager &right );
    
  inline G4bool operator==( const G4EnergyRangeManager &right ) const
    { return ( this == (G4EnergyRangeManager *) &right ); }
    
  inline G4bool operator!=( const G4EnergyRangeManager &right ) const
    { return ( this != (G4EnergyRangeManager *) &right ); }
    
  void RegisterMe( G4HadronicInteraction *a );
    
  G4HadronicInteraction *GetHadronicInteraction(const G4HadProjectile & aHadProjectile, 
                                                G4Nucleus & aTargetNucleus,
						const G4Material *aMaterial,
						const G4Element *anElement ) const;
  // This is the new one to be used.

  G4HadronicInteraction *GetHadronicInteraction(const G4double kineticEnergy,
						const G4Material *aMaterial,
						const G4Element *anElement ) const;
  // This is the old, deprecated one, which will be removed later on.

  std::vector<G4HadronicInteraction*>& GetHadronicInteractionList();
    
  void Dump( G4int verbose = 0 ); 

  void BuildPhysicsTable(const G4ParticleDefinition&);
    
private:
     
  G4int theHadronicInteractionCounter;
  std::vector<G4HadronicInteraction*> theHadronicInteraction;
};

#endif
 
