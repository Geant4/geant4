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
// ClassName:   HadronPhysicsQGSC_BERT
//
// Author: 2007  G.Folger
//             created from HadronPhysicsQGSC originally by J.P. Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------

#ifndef HadronPhysicsQGSC_BERT_h
#define HadronPhysicsQGSC_BERT_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4QGSCPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4QGSCProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4QGSCNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"

class HadronPhysicsQGSC_BERT : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsQGSC_BERT(G4int verbose =1);
    HadronPhysicsQGSC_BERT(const G4String& name,G4bool quasiElastic=true);
    virtual ~HadronPhysicsQGSC_BERT();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4QGSCNeutronBuilder * theQGSCNeutron;
    G4BertiniNeutronBuilder * theBertiniNeutron;
    G4LEPNeutronBuilder * theLEPNeutron;        //needed for capture&fission
    
    G4PiKBuilder * thePiK;
    G4QGSCPiKBuilder * theQGSCPiK;
    G4BertiniPiKBuilder * theBertiniPiK;
    
    G4ProtonBuilder * thePro;
    G4QGSCProtonBuilder * theQGSCPro;    
    G4BertiniProtonBuilder * theBertiniPro;
    
    G4MiscLHEPBuilder * theMiscLHEP;
    
    G4bool QuasiElastic;
};

#endif

