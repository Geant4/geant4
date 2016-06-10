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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsNuBeam
//
// Author: Julia Yarba, FNAL/CD (2014)
// Comment: somewhat "molded" after HadronPhysicsFTFP_BETT
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsNuBeam_h
#define G4HadronPhysicsNuBeam_h 1

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
// specific to NuBeam case
#include "G4QGSPLundStrFragmProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

class G4ComponentGGHadronNucleusXsc;


class G4HadronPhysicsNuBeam : public G4VPhysicsConstructor
{

  public: 
    G4HadronPhysicsNuBeam(G4int verbose =1);
    G4HadronPhysicsNuBeam(const G4String& name, G4bool quasiElastic=false);
    virtual ~G4HadronPhysicsNuBeam();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();

    G4bool QuasiElastic;

    // Simplify handling of TLS data, encapsulate everyhing in a structure
    //
    struct ThreadPrivate { 
    
       G4NeutronBuilder * theNeutrons;
       G4BertiniNeutronBuilder * theBertiniNeutron;
       G4FTFPNeutronBuilder * theFTFPNeutron;
 
       G4PiKBuilder * thePiK;
       G4BertiniPiKBuilder * theBertiniPiK;
       G4FTFPPiKBuilder * theFTFPPiK;
    
       G4ProtonBuilder * thePro;
       G4BertiniProtonBuilder * theBertiniPro;
       G4FTFPProtonBuilder * theFTFPPro;  
       // specific to NuBeam  
       G4QGSPLundStrFragmProtonBuilder * theQGSPPro;
    
       G4HyperonFTFPBuilder * theHyperon;
    
       G4AntiBarionBuilder * theAntiBaryon;
       G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

       G4ComponentGGHadronNucleusXsc * xsKaon;
       G4VCrossSectionDataSet * xsNeutronInelasticXS;
       G4VCrossSectionDataSet * xsNeutronCaptureXS;
    };
    static G4ThreadLocal ThreadPrivate* tpdata;   

};

#endif

