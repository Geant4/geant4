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
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsQGSP_BIC_AllHP_h
#define G4HadronPhysicsQGSP_BIC_AllHP_h 1

#include "G4HadronPhysicsQGSP_BIC.hh"

class G4HadronPhysicsQGSP_BIC_AllHP : public G4HadronPhysicsQGSP_BIC 
{
  public: 
    G4HadronPhysicsQGSP_BIC_AllHP(G4int verbose =1);
    G4HadronPhysicsQGSP_BIC_AllHP(const G4String& name, G4bool quasiElastic=true);
    virtual ~G4HadronPhysicsQGSP_BIC_AllHP() {}

  protected:
    virtual void Neutron() override;
    virtual void Proton() override;
    G4double minBIC_neutron;
    G4double minBIC_proton;
    G4double maxHP_neutron;
    G4double maxHP_proton;
};

#endif

