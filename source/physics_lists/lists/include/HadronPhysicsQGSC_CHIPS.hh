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
// $Id: HadronPhysicsQGSC_CHIPS.hh,v 1.5 2009-04-16 09:26:47 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_CHIPS
//
// Author: 2009  M. Kosov based on HadronPhysicsQGSC_BERT
//
// Modified:
//
//----------------------------------------------------------------------------

#ifndef HadronPhysicsQGSC_CHIPS_h
#define HadronPhysicsQGSC_CHIPS_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
//#include "G4MiscLHEPBuilder.hh"
#include "G4MiscQGSCBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4QGSC_CHIPSPiKBuilder.hh"
//#include "G4BertiniPiKBuilder.hh"

#include "G4QProtonBuilder.hh"
#include "G4QGSC_CHIPSProtonBuilder.hh"
//#include "G4BertiniProtonBuilder.hh"

#include "G4QNeutronBuilder.hh"
#include "G4QGSC_CHIPSNeutronBuilder.hh"
//#include "G4BertiniNeutronBuilder.hh"
//#include "G4LEPNeutronBuilder.hh"

class HadronPhysicsQGSC_CHIPS : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsQGSC_CHIPS(const G4String& name ="hadron",G4bool quasiElastic=true);
    virtual ~HadronPhysicsQGSC_CHIPS();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    G4QNeutronBuilder * theNeutrons;
    G4QGSC_CHIPSNeutronBuilder * theQGSCNeutron;
    //G4BertiniNeutronBuilder * theBertiniNeutron;
    //G4LEPNeutronBuilder * theLEPNeutron;        //needed for capture&fission
    
    G4PiKBuilder * thePiK;
    G4QGSC_CHIPSPiKBuilder * theQGSCPiK;
    //G4BertiniPiKBuilder * theBertiniPiK;
    
    G4QProtonBuilder * thePro;
    G4QGSC_CHIPSProtonBuilder * theQGSCPro;    
    //G4BertiniProtonBuilder * theBertiniPro;
    
    G4MiscQGSCBuilder * theMiscQGSC;
    //G4MiscLHEPBuilder * theMiscLHEP;
    
    G4bool QuasiElastic;
};

#endif

