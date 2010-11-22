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
// $Id: G4EnergyRangeManager.hh,v 1.10 2010-11-22 07:45:33 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Energy Range Manager
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Last modified: 24-Mar-1997
 
#ifndef G4EnergyRangeManager_h
#define G4EnergyRangeManager_h 1
 
#include "G4HadronicInteraction.hh"
 
 class G4EnergyRangeManager 
 {
     
 public:
    
    G4EnergyRangeManager();
 
    ~G4EnergyRangeManager()
    { }
    
    G4EnergyRangeManager(const G4EnergyRangeManager& right);
    
    G4EnergyRangeManager& operator=( const G4EnergyRangeManager &right );
    
 public:
    
    inline G4bool operator==( const G4EnergyRangeManager &right ) const
    { return ( this == (G4EnergyRangeManager *) &right ); }
    
    inline G4bool operator!=( const G4EnergyRangeManager &right ) const
    { return ( this != (G4EnergyRangeManager *) &right ); }
    
    void RegisterMe( G4HadronicInteraction *a );
    
    G4HadronicInteraction *GetHadronicInteraction(
     const G4double kineticEnergy,
     const G4Material *aMaterial,
     const G4Element *anElement ) const;
    
	//private:
    
    inline G4int GetHadronicInteractionCounter() const
    { return theHadronicInteractionCounter; }
    
 private:
     
    enum { MAX_NUMBER_OF_MODELS = 100 };
    
    G4HadronicInteraction *
     theHadronicInteraction[ MAX_NUMBER_OF_MODELS ];
    
    G4int theHadronicInteractionCounter;
    
 };

#endif
 
