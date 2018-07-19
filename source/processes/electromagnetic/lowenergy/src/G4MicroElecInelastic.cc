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
// G4MicroElecInelastic.cc, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//
//          - Inelastic cross-sections of low energy electrons in silicon
//	    for the simulation of heavy ion tracks with theGeant4-DNA toolkit,
//	    NSS Conf. Record 2010, pp. 80-85.
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    NIM B, vol. 288, pp. 66 - 73, 2012.
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for protons and
//	    heavy ions in Si, NIM B, vol. 287, pp. 124 - 129, 2012.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 


#include "G4MicroElecInelastic.hh"
#include "G4DummyModel.hh"
#include "G4SystemOfUnits.hh"

#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MicroElecInelastic::G4MicroElecInelastic(const G4String& processName,
                                           G4ProcessType type):G4VEmProcess (processName, type),
isInitialised(false)
{
  SetProcessSubType(53);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4MicroElecInelastic::~G4MicroElecInelastic()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MicroElecInelastic::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() ||
          &p == G4Proton::Proton()  ||
          &p == G4GenericIon::GenericIonDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MicroElecInelastic::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    G4String name = p->GetParticleName();

    if(!EmModel()) SetEmModel(new G4DummyModel());
    AddEmModel(2, EmModel());   
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MicroElecInelastic::PrintInfo()
{
  // V.I. printout of models is performed by model manager
  //      if this extra printout is needed it should be 
  //      protected by verbosity level
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
