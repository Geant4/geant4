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
//---------------------------------------------------------------------------
//
// ClassName:   G4FTFCaptureAtRest
//
// Author:    Alberto Ribon
//
// Date:      18 October 2011
//
// Modified:
//
// This class provides the nuclear capture at rest process for
// anti-protons, using Fritiof/Precompound model.
//
//----------------------------------------------------------------------------
//

#ifndef G4FTFCaptureAtRest_hh
#define G4FTFCaptureAtRest_hh

#include "globals.hh"
#include "G4ios.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"

class G4TheoFSGenerator;
class G4PreCompoundModel;
class G4GeneratorPrecompoundInterface;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4ExcitationHandler;


class G4FTFCaptureAtRest : public G4VRestProcess {  

public:

  G4FTFCaptureAtRest( const G4String& processName ="FTFNuclearCaptureAtRest" );

  virtual ~G4FTFCaptureAtRest();

  virtual G4bool IsApplicable( const G4ParticleDefinition&  );

  G4VParticleChange* AtRestDoIt( const G4Track& , const G4Step&  ); 

protected:                         

  // Zero mean lifetime
  G4double GetMeanLifeTime( const G4Track& , G4ForceCondition* );

private:

  // Hide assignment operator as private 
  G4FTFCaptureAtRest& operator=( const G4FTFCaptureAtRest &right );

  // Copy constructor
  G4FTFCaptureAtRest( const G4FTFCaptureAtRest& );

  // Useful method to print out information in case of problems
  void DumpState( const G4Track& , const G4String& );
 
  G4TheoFSGenerator * theModel;
  G4PreCompoundModel * thePreEquilib;
  G4GeneratorPrecompoundInterface * theCascade;
  G4FTFModel * theStringModel;
  G4ExcitedStringDecay * theStringDecay;
  G4LundStringFragmentation * theLund;
  G4ExcitationHandler * theHandler;

  G4double theMin;
  G4double theMax;

};


inline G4double G4FTFCaptureAtRest::
GetMeanLifeTime( const G4Track& , G4ForceCondition* ) {
  return 0.0;
}


#endif

