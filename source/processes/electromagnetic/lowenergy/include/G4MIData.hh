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
// gpaterno, March 2019
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4MIData_h
#define G4MIData_h 1

#include "globals.hh"
#include "G4VMaterialExtension.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MIData : public G4VMaterialExtension {

public:  
  explicit G4MIData(const G4String&);
  virtual ~G4MIData();
    
public:
  void Print() const override 
  {G4cout << "Molecular Interference data for Rayleigh scattering" << G4endl;};
        
  void SetFilenameFF(const G4String& filenameff) {fFilenameFF = filenameff;};
  void SetFilenameCS(const G4String& filenamecs) {fFilenameCS = filenamecs;};
  void SetMolWeight(const G4double mw) {fMolWeight = mw;};
  
  const G4String& GetFilenameFF() {return fFilenameFF;};
  const G4String& GetFilenameCS() {return fFilenameCS;};
  const G4double& GetMolWeight() {return fMolWeight;};
  
private:
  G4String fFilenameFF;
  G4String fFilenameCS;
  G4double fMolWeight;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif







