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
#ifndef F01RunAction_h
#define F01RunAction_h

#include "G4UserRunAction.hh"

class G4ParticleDefinition;
class G4Transportation;
class G4CoupledTransportation;

class G4Run;

class F01RunAction: public G4UserRunAction {

public:

  F01RunAction();
  virtual ~F01RunAction();
      
  virtual void BeginOfRunAction( const G4Run* aRun );    
  virtual void EndOfRunAction( const G4Run* aRun );    

  // Helper method to change the Transportation's 'looper' parameters 
  void ChangeLooperParameters(const G4ParticleDefinition* particleDef );

  // Helper method to find the Transportation process for a particle type
  G4Transportation*
     FindTransportation( const G4ParticleDefinition * particleDef,
                         bool reportError= true );

public:
  void     SetNumberOfTrials( G4int val )   { fNumberOfTrials  =  val; }
  void     SetWarningEnergy( G4double val )   { fWarningEnergy = val; }
  void     SetImportantEnergy( G4double val ) { fImportantEnergy = val; }
  G4int    GetNumberOfTrials() { return fNumberOfTrials; }
  G4double GetWarningEnergy()  { return fWarningEnergy; }
  G4double GetImportantEnergy() { return fImportantEnergy; }
   
private:

  // Values for initialising 'loopers' parameters of Transport process
  G4int    fNumberOfTrials  =  0;    // Default will not overwrite
  G4double fWarningEnergy   = -1.0;  // Default values - non operational
  G4double fImportantEnergy = -1.0;  // Default - will not overwrite

  G4int    fVerboseLevel = 0;
};

#endif
