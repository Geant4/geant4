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
/// \file electromagnetic/TestEm11/include/Run.hh
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
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);  

    void AddEdep (G4int i, G4double e);
    void AddTotEdep     (G4double e);
    void AddTrackLength (G4double t);
    void AddProjRange   (G4double x);
    void AddStepSize    (G4int nb, G4double st);
    void AddTrackStatus (G4int i);
    
    void SetCsdaRange (G4int i, G4double value);
    void SetXfrontNorm(G4int i, G4double value);
                                      
    G4double GetCsdaRange (G4int i);
    G4double GetXfrontNorm(G4int i);   
            
    virtual void Merge(const G4Run*);
    void EndOfRun();
    
  private:
    DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin; 

    G4double   fTrackLen,  fTrackLen2;
    G4double   fProjRange, fProjRange2;
    G4int      fNbOfSteps, fNbOfSteps2;
    G4double   fStepSize,  fStepSize2;
    G4int      fStatus[3];
    
    G4double   fEdeposit[kMaxAbsor];
    G4double   fEmin[kMaxAbsor], fEmax[kMaxAbsor];
    G4double   fTotEdep[3];
    G4double   fCsdaRange[kMaxAbsor];
    G4double   fXfrontNorm[kMaxAbsor];    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

