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
// ClassName:   
//
// Author: 7 Nov 2017 Tatsumi Koi
//   created from G4HadronPhysicsShielding
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsShieldingLEND_h
#define G4HadronPhysicsShieldingLEND_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VPhysicsConstructor.hh"


class G4HadronPhysicsShieldingLEND : public G4VPhysicsConstructor
{
  public:
    explicit G4HadronPhysicsShieldingLEND(G4int verbose=1);
    explicit G4HadronPhysicsShieldingLEND(const G4String& name, G4bool );
    explicit G4HadronPhysicsShieldingLEND(const G4String& name, G4int verbose=1,
                                          G4double minFTFPEnergy=9.5*CLHEP::GeV, G4double maxBertiniEnergy=9.9*CLHEP::GeV);
    virtual ~G4HadronPhysicsShieldingLEND();

  public: 
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    void UseLEND( G4String ss="" ){ useLEND_=true; evaluation_=ss; };
    void UnuseLEND(){ useLEND_=false; };

  private:
    virtual void CreateModels();
    virtual void Neutron();
    virtual void Proton();
    virtual void Pion();
    virtual void Kaon();
    virtual void Others();
    virtual void DumpBanner();
    //This contains extra configurataion specific to this PL
    virtual void ExtraConfiguration();

    G4bool useLEND_;
    G4String evaluation_;

    const G4double minFTFPEnergy_;
    const G4double maxBertiniEnergy_;
    const G4double minNonHPNeutronEnergy_;
};

#endif

