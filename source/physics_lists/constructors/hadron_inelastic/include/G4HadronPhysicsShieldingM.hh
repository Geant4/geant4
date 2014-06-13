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
// $Id: G4HadronPhysicsShieldingM.hh 73040 2013-08-15 09:36:57Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2014 K. Genser based on G4HadronPhysicsShielding
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsShieldingM_h
#define G4HadronPhysicsShieldingM_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4FTFPProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4NeutronHPBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"


class G4HadronPhysicsShieldingM : public G4VPhysicsConstructor
{
  public: 
    //G4HadronPhysicsShieldingM(G4int verbose =1,G4bool blend=false);
    explicit G4HadronPhysicsShieldingM(G4int verbose=1);
    explicit G4HadronPhysicsShieldingM(const G4String& name);
    explicit G4HadronPhysicsShieldingM(const G4String& name, G4int verbose=1,
                                       G4double minFTFPEnergy=9.5*GeV, G4double maxBertiniEnergy=9.9*GeV);
    virtual ~G4HadronPhysicsShieldingM();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    void UseLEND( G4String ss="" ){useLEND_=true;evaluation_=ss;};
    void UnuseLEND(){useLEND_=false;};

  private:
    void CreateModels();
    
    struct ThreadPrivate { 
      G4NeutronBuilder * theNeutrons;
      //G4NeutronHPBuilder * theHPNeutron;
      G4VNeutronBuilder * theLENeutron;
      G4BertiniNeutronBuilder * theBertiniNeutron;
      G4FTFPNeutronBuilder * theFTFPNeutron;
 
      G4PiKBuilder * thePiK;
      G4BertiniPiKBuilder * theBertiniPiK;
      G4FTFPPiKBuilder * theFTFPPiK;
    
      G4ProtonBuilder * thePro;
      G4BertiniProtonBuilder * theBertiniPro;
      G4FTFPProtonBuilder * theFTFPPro;    

      G4HyperonFTFPBuilder * theHyperon;
    
      G4AntiBarionBuilder * theAntiBaryon;
      G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

      G4VCrossSectionDataSet * theChipsKaonMinus;
      G4VCrossSectionDataSet * theChipsKaonPlus;
      G4VCrossSectionDataSet * theChipsKaonZero;
 
      G4VCrossSectionDataSet * theBGGxsNeutron;
      G4VCrossSectionDataSet * theNeutronHPJENDLHEInelastic;
      G4VCrossSectionDataSet * theBGGxsProton;

      G4VCrossSectionDataSet * xsNeutronCaptureXS;
    };
    static G4ThreadLocal ThreadPrivate* tpdata;

    // G4bool QuasiElastic;
    G4bool useLEND_;
    G4String evaluation_;

    const G4double minFTFPEnergy_;
    const G4double maxBertiniEnergy_;
    const G4double minNonHPNeutronEnergy_;

};

#endif
