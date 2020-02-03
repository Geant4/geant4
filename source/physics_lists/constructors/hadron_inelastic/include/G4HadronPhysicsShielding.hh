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
// Author: 2007  tatsumi Koi, Gunter Folger
//   created from G4HadronPhysicsFTFP_BERT
// Modified:
// 2019.08.01 A.Ribon replaced explicit numbers for the energy transition
//                    region with values taken from G4HadronicParameters
// 2014.08.05 K.L.Genser added provisions for modifing the Bertini to
//            FTF transition energy region
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsShielding_h
#define G4HadronPhysicsShielding_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VPhysicsConstructor.hh"
#include "G4HadronicParameters.hh"


class G4HadronPhysicsShielding : public G4VPhysicsConstructor
{
  public: 
    explicit G4HadronPhysicsShielding(G4int verbose=1);
    explicit G4HadronPhysicsShielding(const G4String& name, G4bool );
    explicit G4HadronPhysicsShielding(const G4String& name, G4int verbose=1,
      G4double minFTFPEnergy=G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade(),
      G4double maxBertiniEnergy=G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade());

    virtual ~G4HadronPhysicsShielding();

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

