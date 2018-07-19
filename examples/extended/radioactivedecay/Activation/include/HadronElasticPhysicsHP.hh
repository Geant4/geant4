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
/// \file HadronElasticPhysicsHP.hh
/// \brief Definition of the HadronElasticPhysicsHP class
//
// $Id: HadronElasticPhysicsHP.hh 71037 2013-06-10 09:20:54Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HadronElasticPhysicsHP_h
#define HadronElasticPhysicsHP_h 1

#include "globals.hh"
#include "G4HadronElasticPhysics.hh"

class NeutronHPMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HadronElasticPhysicsHP : public G4HadronElasticPhysics
{
  public: 
    HadronElasticPhysicsHP(G4int ver = 1); 
   ~HadronElasticPhysicsHP();

  public: 
    virtual void ConstructProcess();
    
  public:
    void SetThermalPhysics(G4bool flag) {fThermal = flag;};
      
  private:
    G4bool                  fThermal;
    NeutronHPMessenger*     fNeutronMessenger;          
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

