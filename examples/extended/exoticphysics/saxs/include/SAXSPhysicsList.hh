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
/// \file SAXSPhysicsList.hh
/// \brief Implementation of the SAXSPhysicsList class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SAXSPhysicsList_h
#define SAXSPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include <vector>

class G4VPhysicsConstructor;
class SAXSPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Physics list.
/// It includes various EM constructors, which can be selected through macro.
/// By defualt "standard" Penelope physics is used. In order to activate 
/// Molecular interference (MI) effects in Rayleigh scattering, choose 
/// G4EmPenelopeMI PhysicsList with fUseMIFlag variable set as true (default).

class SAXSPhysicsList : public G4VUserPhysicsList
{
public:
  SAXSPhysicsList();
  virtual ~SAXSPhysicsList();

  void ConstructParticle() override;
  void ConstructProcess() override;
   
  //for the Messenger 
  void SetDefaultCutsValue(G4double);
  void SelectPhysicsList(const G4String& name);
  void SetUseMIFlag(G4bool val){fUseMIFlag = val;};
  G4bool GetUseMIFlag(){return fUseMIFlag;};
 
  SAXSPhysicsList & operator = (const SAXSPhysicsList &right) = delete;
  SAXSPhysicsList(const SAXSPhysicsList&) = delete;

private:
  
  G4VPhysicsConstructor* fParticleList; 
  G4VPhysicsConstructor* fEmPhysicsList;    
  
  G4bool fUseMIFlag; 
      
  SAXSPhysicsListMessenger* fPMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif 

