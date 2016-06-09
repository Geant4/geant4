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
// $Id: G4EmConfigurator.hh,v 1.2 2008/11/21 12:30:29 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4EmConfigurator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 14.07.2008
//
// Modifications:
//
// Class Description:
//
// This class provides configuration EM models for 
// particles/processes/regions
//

// -------------------------------------------------------------------
//

#ifndef G4EmConfigurator_h
#define G4EmConfigurator_h 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmConfigurator 
{
public: 
  
  G4EmConfigurator();
 
  ~G4EmConfigurator();

  // Add EM model to the list of extra models potentially to be 
  // declared for the G4Region and energy interval
  // 
  void AddExtraEmModel(const G4String& particleName,
		       G4VEmModel*, G4VEmFluctuationModel* fm = 0); 

  // Declare EM model for particle type and process to 
  // be active for the G4Region and energy interval
  // The model should be previously added to the configurator 
  // or be "dummy"
  // 
  void AddModelForRegion(const G4String& particleName,
                         const G4String& processName,
                         const G4String& modelName,
                         const G4String& regionName = "",
                         G4double emin = 0.0,
                         G4double emax = DBL_MAX,
                         const G4String& flucModelName = "");

  // Set EM model for particle type and process to 
  // be active for the G4Region and energy interval
  //
  void SetExtraEmModel(const G4String& particleName,
		       const G4String& processName,
		       G4VEmModel*,
		       const G4String& regionName = "",
		       G4double emin = 0.0,
		       G4double emax = DBL_MAX,
		       G4VEmFluctuationModel* fm = 0); 

  // Add all previously declared models to corresponding processes
  void AddModels();

private:

  void SetModelForRegion(const G4String& particleName,
                         const G4String& processName,
                         const G4String& modelName,
                         const G4String& regionName,
                         const G4String& flucModelName,
                         G4double emin,
                         G4double emax);

  // hide assignment operator
  G4EmConfigurator & operator=(const G4EmConfigurator &right);
  G4EmConfigurator(const G4EmConfigurator&);

  std::vector<G4String> particles;  
  std::vector<G4String> processes;  
  std::vector<G4String> models;  
  std::vector<G4String> regions;  
  std::vector<G4String> flucModels;  
  std::vector<G4double> lowEnergy;
  std::vector<G4double> highEnergy;
  
  std::vector<G4String> particleList;  
  std::vector<G4VEmModel*> modelList;
  std::vector<G4VEmFluctuationModel*> flucModelList;

  G4int index;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif








