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
// $Id: HadronPhysicsFTFP_BERT.hh,v 1.5 2010-06-03 10:42:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007  Gunter Folger
//   created from HadronPhysicsFTFP
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 19.06.2008 G.Folger: change default for QE to NOT use Chips QE
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsFTFP_BERT_h
#define HadronPhysicsFTFP_BERT_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscCHIPSBuilder.hh"

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
#include "G4LEPNeutronBuilder.hh"

class HadronPhysicsFTFP_BERT : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsFTFP_BERT(G4int verbose =1);
    HadronPhysicsFTFP_BERT(const G4String& name, G4bool quasiElastic=false);
    virtual ~HadronPhysicsFTFP_BERT();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4HadronicProcess* FindInelasticProcess(const G4ParticleDefinition*);
    
    G4NeutronBuilder * theNeutrons;
    G4BertiniNeutronBuilder * theBertiniNeutron;
    G4FTFPNeutronBuilder * theFTFPNeutron;
    G4LEPNeutronBuilder * theLEPNeutron;        //needed for capture&fission
 
    G4PiKBuilder * thePiK;
    G4BertiniPiKBuilder * theBertiniPiK;
    G4FTFPPiKBuilder * theFTFPPiK;
    
    G4ProtonBuilder * thePro;
    G4BertiniProtonBuilder * theBertiniPro;
    G4FTFPProtonBuilder * theFTFPPro;    
    
    G4MiscCHIPSBuilder * theMiscCHIPS;
    
    G4bool QuasiElastic;
    G4VCrossSectionDataSet * theCHIPSInelastic;
};

#endif

