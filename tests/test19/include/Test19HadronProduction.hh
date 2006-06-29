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
// ---------------------------------------------------------------------
// Short description of the G4 Class made by M.K. (asFarAsHeNnderstands)
// =====================================================================
// >>> StandardG4Process wrepper for the Hadronic Model Collection <<<
// This class is necessary only for the navigation through all models
// of the Hadronic collection AND for conversion of the special classes
// of the Hadronic Result (list of the secondary particles) to the
// standard G4VParticleChange form. One should not use this class
// for the simple ElectroWeak Processes, which are using the only
// event Generator CHIPS and produce the Result in the
// standard G4VParticleChange form. -- 29.01.2004. -- M.K.
//
// ***Important*** It is working only for the G4VDiscreteProcess
//                 (not for G4VRestProcess - @@ Create similar class)
//
// ======================================================================================
//
//      ---------- Test19HadronProduction -------
//    Converted from Test19 to Test19 by Mikhail Kossov, 29 Jan 2005 
//
//---------------------------------------------------------------------------------------

#ifndef Test19HadronProduction_Test19HadronProduction_h
#define Test19HadronProduction_Test19HadronProduction_h 1

#include "Test19VSecondaryGenerator.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4HadFinalState.hh"

class G4Step;
class G4Track;
class G4VParticleChange;

class Test19HadronProduction : public G4VDiscreteProcess
{
public:

  Test19HadronProduction(const G4String& processName = "HadronProduction" );
  ~Test19HadronProduction();

  void SetSecondaryGenerator(Test19VSecondaryGenerator*);
  G4double PostStepGetPhysicalInteractionLength(const             G4Track& track,
                                                G4double          previousStepSize,
                                                G4ForceCondition* condition);
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
  // @@This is a bad stile to declare and to implement in one line -> make inline
  G4bool IsApplicable(const G4ParticleDefinition&) {return true;};
  G4double GetMass() {return theGenerator->GetMass();};

//protected:
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) {return DBL_MAX;};

private:
  // hide assignment operator as private
  Test19HadronProduction(const Test19HadronProduction&);
  Test19HadronProduction& operator = (const Test19HadronProduction &right);

  // Navigation tool for the Model Collection of the Hadronic package
  void InitializeMe();

  // Body
  Test19VSecondaryGenerator* theGenerator; // The selected model from the Collection
  G4VParticleChange          theChange;    // Standard G4VProcess output Prototype
};

#endif  // Test19HadronProduction_Test19HadronProduction_h

