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
#ifndef G4HadronPhysicsQGSP_BIC_AllHP_h
#define G4HadronPhysicsQGSP_BIC_AllHP_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4ProtonPHPBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4NeutronPHPBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

#include "G4HadronicParameters.hh"


class G4ComponentGGHadronNucleusXsc;


class G4HadronPhysicsQGSP_BIC_AllHP : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsQGSP_BIC_AllHP(G4int verbose =1);
    G4HadronPhysicsQGSP_BIC_AllHP(const G4String& name, G4bool quasiElastic=true);
    virtual ~G4HadronPhysicsQGSP_BIC_AllHP();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();

    const G4double maxFTFP = G4HadronicParameters::Instance()->GetMaxEnergyTransitionQGS_FTF();
    const G4double minFTFP = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
    const G4double maxBIC  = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
    const G4double maxBERT = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
    const G4double minBIC_neutron = 19.9*CLHEP::MeV;
    const G4double maxHP_neutron =  20.0*CLHEP::MeV;
    const G4double minBIC_proton = 199.0*CLHEP::MeV;
    const G4double maxHP_proton =  200.0*CLHEP::MeV;

    struct ThreadPrivate {
      G4NeutronBuilder * theNeutrons;
      G4FTFPNeutronBuilder * theFTFPNeutron;
      G4QGSPNeutronBuilder * theQGSPNeutron;
      G4BinaryNeutronBuilder * theBinaryNeutron;
      G4NeutronPHPBuilder * thePHPNeutron;
    
      G4PiKBuilder * thePiK;
      G4FTFPPiKBuilder * theFTFPPiK;
      G4QGSPPiKBuilder * theQGSPPiK;
      G4BertiniPiKBuilder * theBertiniPiK;
    
      G4ProtonBuilder * thePro;
      G4FTFPProtonBuilder * theFTFPPro;
      G4QGSPProtonBuilder * theQGSPPro;
      G4BinaryProtonBuilder * theBinaryPro;
      G4ProtonPHPBuilder * thePHPProton;

      G4HyperonFTFPBuilder * theHyperon;

      G4AntiBarionBuilder * theAntiBaryon;
      G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;
    };
    static G4ThreadLocal ThreadPrivate* tpdata;

};

#endif
