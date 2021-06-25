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
/// \file medical/dna/range/include/Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4Event.hh"

#include "G4THitsMap.hh"
#include <vector>

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(const DetectorConstruction* /*detector*/);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);  

    G4double GetEdep()          const {return fEdeposit;};

    virtual void RecordEvent(const G4Event*);
    G4THitsMap<G4double>* GetHitsMap(){return fRunMap;}
    G4THitsMap<G4double>* GetHitsMap(const G4String& detName, 
                                     const G4String& colName);
    G4THitsMap<G4double>* GetHitsMap(const G4String& fullName);
    void DumpAllScorer();
        
    virtual void Merge(const G4Run*);
    void EndOfRun();
    
  private:
    //const DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin; 
    G4double  fEdeposit; 
       
    G4String fCollName;
    //G4int fCollID;
    G4THitsMap<G4double>* fRunMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

