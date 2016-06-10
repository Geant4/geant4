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
// $Id: G4HadronPhysicsQGS_BIC.hh 76703 2013-11-14 10:29:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGS_BIC
//
// Author: 2007 Gunter Folger
//     created from G4HadronPhysicsQGSP_BIC  by H.P.Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsQGS_BIC_h
#define G4HadronPhysicsQGS_BIC_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PionBuilder.hh"
#include "G4BinaryPionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4FTFBinaryPionBuilder.hh"
#include "G4QGSBinaryPionBuilder.hh"

#include "G4KaonBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFBinaryKaonBuilder.hh"
#include "G4QGSBinaryKaonBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4FTFBinaryProtonBuilder.hh"
#include "G4QGSBinaryProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFBinaryNeutronBuilder.hh"
#include "G4QGSBinaryNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"


class G4HadronPhysicsQGS_BIC : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsQGS_BIC(G4int verbose =1);
    G4HadronPhysicsQGS_BIC(const G4String& name, G4bool quasiElastic=true);
    virtual ~G4HadronPhysicsQGS_BIC();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();

    struct ThreadPrivate {
      G4NeutronBuilder * theNeutrons;
      G4FTFBinaryNeutronBuilder * theFTFBinaryNeutron;
      G4QGSBinaryNeutronBuilder * theQGSBinaryNeutron;
      G4BinaryNeutronBuilder * theBinaryNeutron;

      G4PionBuilder * thePion;
      G4BinaryPionBuilder * theBinaryPion;
      G4BertiniPionBuilder * theBertiniPion;
      G4FTFBinaryPionBuilder * theFTFBinaryPion;
      G4QGSBinaryPionBuilder * theQGSBinaryPion;

      G4KaonBuilder * theKaon;
      G4BertiniKaonBuilder * theBertiniKaon;
      G4FTFBinaryKaonBuilder * theFTFBinaryKaon;
      G4QGSBinaryKaonBuilder * theQGSBinaryKaon;

      G4ProtonBuilder * thePro;
      G4FTFBinaryProtonBuilder * theFTFBinaryPro;
      G4QGSBinaryProtonBuilder * theQGSBinaryPro;
      G4BinaryProtonBuilder * theBinaryPro;

      G4HyperonFTFPBuilder * theHyperon;

      G4AntiBarionBuilder * theAntiBaryon;
      G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

      G4VCrossSectionDataSet * xsNeutronInelasticXS;
      G4VCrossSectionDataSet * xsNeutronCaptureXS;
    };
    static G4ThreadLocal ThreadPrivate* tpdata;

    // G4bool QuasiElastic;
};

#endif

