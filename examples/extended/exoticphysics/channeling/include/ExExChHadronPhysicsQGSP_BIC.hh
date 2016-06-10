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
// $Id: ExExChHadronPhysicsQGSP_BIC.hh 76703 2013-11-14 10:29:11Z gcosmo $
//

#ifndef ExExChHadronPhysicsQGSP_BIC_h
#define ExExChHadronPhysicsQGSP_BIC_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4FTFPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

#include "ExExChProtonBuilder.hh"
#include "ExExChPiKBuilder.hh"
#include "ExExChHyperonFTFPBuilder.hh"
#include "ExExChAntiBarionBuilder.hh"


class ExExChHadronPhysicsQGSP_BIC : public G4VPhysicsConstructor
{
  public: 
    ExExChHadronPhysicsQGSP_BIC(G4int verbose =1);
    ExExChHadronPhysicsQGSP_BIC(const G4String& name,
                                      G4bool quasiElastic=true);
    virtual ~ExExChHadronPhysicsQGSP_BIC();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();

    struct ThreadPrivate {
      G4NeutronBuilder * theNeutrons;
      G4FTFPNeutronBuilder * theFTFPNeutron;
      G4QGSPNeutronBuilder * theQGSPNeutron;
      G4BinaryNeutronBuilder * theBinaryNeutron;
    
      ExExChPiKBuilder * thePiK;
      G4FTFPPiKBuilder * theFTFPPiK;
      G4QGSPPiKBuilder * theQGSPPiK;
      G4BertiniPiKBuilder * theBertiniPiK;
    
      ExExChProtonBuilder * thePro;
      G4FTFPProtonBuilder * theFTFPPro;
      G4QGSPProtonBuilder * theQGSPPro; 
      G4BinaryProtonBuilder * theBinaryPro;

      ExExChHyperonFTFPBuilder * theHyperon;

      ExExChAntiBarionBuilder * theAntiBaryon;
      G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

      G4VCrossSectionDataSet * xsNeutronInelasticXS;
      G4VCrossSectionDataSet * xsNeutronCaptureXS;
    };
    static G4ThreadLocal ThreadPrivate* tpdata;

    // G4bool QuasiElastic;
};

#endif

