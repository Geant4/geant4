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
// ClassName:   HadronPhysicsQGSP_BERT_NOLEP
//
// Author: 2008 G.Folger
//       created from HadronPhysicsQGSP_BERT
// Modified:
// 15.10.2005 G.Folger: first version
// 31.10.2012 A.Ribon: Use G4MiscBuilder
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsQGSP_BERT_NOLEP_h
#define HadronPhysicsQGSP_BERT_NOLEP_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"

class HadronPhysicsQGSP_BERT_NOLEP : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsQGSP_BERT_NOLEP(G4int verbose =1);
     HadronPhysicsQGSP_BERT_NOLEP(const G4String& name, G4bool quasiElastic=true);
   virtual ~HadronPhysicsQGSP_BERT_NOLEP();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void SetQuasiElastic(G4bool value) {QuasiElastic = value;}; 
    void SetProjectileDiffraction(G4bool value) {ProjectileDiffraction = value;}; 

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4LEPNeutronBuilder * theLEPNeutron;
    G4QGSPNeutronBuilder * theQGSPNeutron;
    G4BertiniNeutronBuilder * theBertiniNeutron;
    
    G4PiKBuilder * thePiK;
    G4QGSPPiKBuilder * theQGSPPiK;
    G4BertiniPiKBuilder * theBertiniPiK;
    
    G4ProtonBuilder * thePro;
    G4QGSPProtonBuilder * theQGSPPro; 
    G4BertiniProtonBuilder * theBertiniPro;
    
    G4MiscBuilder * theMisc;
    
    G4bool QuasiElastic;
    G4bool ProjectileDiffraction;
};

#endif

