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
/// \file electromagnetic/TestEm12/include/Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction* detector);
   ~Run() override = default;

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);  

    void AddEdep (G4double e);
    void AddTrackLength (G4double t);
    void AddProjRange   (G4double x);
    void AddStepSize    (G4int nb, G4double st);
    
    void     SetCsdaRange (G4double value);                                 
    G4double GetCsdaRange();
            
    void Merge(const G4Run*) override;
    void EndOfRun();
    
  private:
    DetectorConstruction*  fDetector = nullptr;
    G4ParticleDefinition*  fParticle = nullptr;
    G4double  fEkin = 0.; 
       
    G4double   fEdeposit  = 0., fEdeposit2  = 0.;
    G4double   fTrackLen  = 0., fTrackLen2  = 0.;
    G4double   fProjRange = 0., fProjRange2 = 0.;
    G4int      fNbOfSteps = 0,  fNbOfSteps2 = 0;
    G4double   fStepSize  = 0., fStepSize2  = 0.;
    
    G4double   fCsdaRange = 0.;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

