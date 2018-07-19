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
// $Id: G4HadronPhysicsShielding.hh 101813 2016-11-30 18:05:43Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007  tatsumi Koi, Gunter Folger
//   created from G4HadronPhysicsFTFP_BERT
// Modified:
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
#include "G4NeutronPHPBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

class G4ComponentGGHadronNucleusXsc;


class G4HadronPhysicsShielding : public G4VPhysicsConstructor
{
  public: 
    //G4HadronPhysicsShielding(G4int verbose =1,G4bool blend=false);
    explicit G4HadronPhysicsShielding(G4int verbose=1);
    explicit G4HadronPhysicsShielding(const G4String& name, G4bool );
    explicit G4HadronPhysicsShielding(const G4String& name, G4int verbose=1,
                                      G4double minFTFPEnergy=9.5*CLHEP::GeV, G4double maxBertiniEnergy=9.9*CLHEP::GeV);
    virtual ~G4HadronPhysicsShielding();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    void UseLEND( G4String ss="" ){useLEND_=true;evaluation_=ss;};
    void UnuseLEND(){useLEND_=false;};

  private:
    void CreateModels();
    
    struct ThreadPrivate { 
      G4NeutronBuilder * theNeutrons;
      //G4NeutronPHPBuilder * theHPNeutron;
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

      G4ComponentGGHadronNucleusXsc * xsKaon;
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

