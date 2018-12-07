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
/// \file medical/electronScattering/include/Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4ParticleDefinition.hh"

class DetectorConstruction;
class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle);
    void EndOfRun();
    void SumFluence(G4double, G4double);
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);
    virtual void Merge(const G4Run*);



  private:
    void InitFluence ();
    void ComputeFluenceError();
    void PrintFluence(G4int);

   DetectorConstruction*  fDetector;
   G4ParticleDefinition* fParticle;
   G4double              fEnergy;


    //for fluence computation

   G4int                   fNbBins;
   G4double                fDr;
   std::vector<G4double>   fluence;
   std::vector<G4double>   fluence1;        //normalized fluence    
   std::vector<G4double>   fluence2;        //rms on norm. fl
   std::vector<G4int>      fNbEntries;      //entries per bin       

};
#endif
