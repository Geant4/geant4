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
#ifndef G4DiscreteScatteringProcess_HH
#define G4DiscreteScatteringProcess_HH

#include "G4VEmProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DiscreteScatteringProcess : public G4VEmProcess{
        
public:
        
  G4DiscreteScatteringProcess(G4int iNumAngles=1);
  
  virtual ~G4DiscreteScatteringProcess();
  
  virtual void InitialiseProcess(const G4ParticleDefinition*);
  
  virtual bool IsApplicable(const G4ParticleDefinition& p);
  // This will need to check whether the particle is something 
  // our particle works for. Right now, protons.
  
  virtual void PrintInfo();
  
private:
     
  bool  fIsInitialised;    
  G4int fNumAngles;
};




#endif


