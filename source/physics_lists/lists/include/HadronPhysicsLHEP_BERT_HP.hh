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
// $Id: HadronPhysicsLHEP_BERT_HP.hh,v 1.1 2006/10/31 11:35:03 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:  HadronPhysicsLHEP_BERT_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsLHEP_BERT_HP_h
#define HadronPhysicsLHEP_BERT_HP_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4LHEPPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LHEPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LHEPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4NeutronHPBuilder.hh"

class HadronPhysicsLHEP_BERT_HP : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsLHEP_BERT_HP(const G4String& name ="hadron");
    virtual ~HadronPhysicsLHEP_BERT_HP();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4NeutronBuilder * theNeutrons;
    G4LHEPNeutronBuilder * theLHEPNeutron;
    G4BertiniNeutronBuilder * theBertiniNeutron;
    G4NeutronHPBuilder * theHPNeutron;
    
    G4PiKBuilder * thePiK;
    G4LHEPPiKBuilder * theLHEPPiK;
    
    G4ProtonBuilder * thePro;
    G4LHEPProtonBuilder * theLHEPPro;
    G4BertiniProtonBuilder * theBertiniPro;
    
    G4MiscLHEPBuilder * theMiscLHEP;
};

#endif

