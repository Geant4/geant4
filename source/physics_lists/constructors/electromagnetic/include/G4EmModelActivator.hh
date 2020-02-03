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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmModelActivator
//
// Author:      V.Ivanchenko 01.06.2015
//
// Organisation:   G4AI
// Contract:       ESA contract 4000107387/12/NL/AT
// Customer:       ESA/ESTEC
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4EmModelActivator_h
#define G4EmModelActivator_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmParameters;
class G4ProcessManager;
class G4ParticleDefinition;
class G4EmConfigurator;
class G4VMscModel;

class G4EmModelActivator 
{
public:

  explicit G4EmModelActivator(const G4String& emphys="");

private:

  void ActivatePAI();

  void ActivateMicroElec();

  void ActivateEmOptions();

  void FindOrAddProcess(const G4ParticleDefinition*, const G4String&);

  void AddStandardScattering(const G4ParticleDefinition*, G4EmConfigurator*,
                             G4VMscModel*, const G4String&, 
                             G4double, G4double);

  G4bool HasMsc(G4ProcessManager*) const;

  G4VMscModel* GetGSModel();

  G4EmModelActivator & operator=(const G4EmModelActivator &right) = delete;
  G4EmModelActivator(const G4EmModelActivator&) = delete;

  G4EmParameters* theParameters;
  G4String baseName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

