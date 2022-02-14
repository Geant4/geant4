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
/// \file 
/// \brief Definition of the G4SteppingVerboseWithUnits class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Stepping Verbose with units for all the applicable double values
//  This class is ported from TestEm2 extended example

//  Original author : Michel Maire (LAPP)
//  Porting with addition of UI command : Makoto Asai (SLAC) Feb.23.2021
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4SteppingVerboseWithUnits_hh
#define G4SteppingVerboseWithUnits_hh 1

#include "G4SteppingVerbose.hh"

class G4GenericMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4SteppingVerboseWithUnits : public G4SteppingVerbose
{
  public:   

    G4SteppingVerboseWithUnits(G4int precision = 4);
   ~G4SteppingVerboseWithUnits() override;
 
    G4VSteppingVerbose* Clone() override
      { return new G4SteppingVerboseWithUnits; }

    void SetManager(G4SteppingManager* const) override;
   
    void TrackingStarted() override;
    void StepInfo() override;

    void AtRestDoItInvoked() override;
    void AlongStepDoItAllDone()override;
    void PostStepDoItAllDone()override;
    void AlongStepDoItOneByOne()override;
    void PostStepDoItOneByOne()override;
    void DPSLStarted() override;
    void DPSLUserLimit() override;
    void DPSLPostStep() override;
    void DPSLAlongStep() override;
    void VerboseTrack() override;
    void VerboseParticleChange() override;
    void ShowStep() const override;
  
  private:

    G4int fprec;  
    G4GenericMessenger* fmessenger = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
